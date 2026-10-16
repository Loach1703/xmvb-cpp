#include "vbscf/structures/assembly/action.hpp"

#include <algorithm>
#include <numeric>

#include "parallel/openmp.hpp"
#include "core/timing.hpp"
#include "mcscf/ci/product_support.hpp"
#include "mcscf/ci/sigma.hpp"
#include "vbscf/structures/assembly/column_tiles.hpp"
#include "vbscf/structures/orthogonal_ci/exterior_transform.hpp"

namespace xmvb::vb {

StructureResponseSweep StructureAction::product_jet(const StructureAdjointState& state,
    const Eigen::MatrixXd* response,const Eigen::MatrixXd* delta_upper,
    const mcscf::DirectCiSigmaAction* delta_sigma,const ExteriorOrbitalTransform& alpha,
    const ExteriorOrbitalTransform& beta,const mcscf::DirectCiSigmaAction& sigma,
    bool collect_images,StructureIntegralDirection* diagnostic) const {
  const bool tangent=response || delta_upper || delta_sigma;
  const int width=tangent ? 2:1,n=sigma.n_orbitals(),d=sigma.channel_kernel().rows();
  const int states=state.structure_coefficients.cols();
  StructureResponseSweep result{
      {Eigen::MatrixXd::Zero(collect_images ? n_structures_:0,states),
       Eigen::MatrixXd::Zero(collect_images ? n_structures_:0,states)},
      Eigen::MatrixXd::Zero(n,n),Eigen::MatrixXd::Zero(d,d),Eigen::MatrixXd::Zero(n,n)};
  // Residency only controls batching. It never selects a different wave function.
  const std::size_t cache_bytes=state.product_images ? state.product_images->budget:
      detail::kStructureTileWorkspaceBytes/4;
  const int extent=std::max<std::size_t>(1,std::min<std::size_t>(n_unique_beta_,
      std::min<std::size_t>(8*cooperative_openmp_thread_count(),
          // Three output components plus compact row indices must fit even
          // for highly unequal spin spaces; the bound is shared by accepted
          // and directional calls so diagnostic capture does not change it.
          cache_bytes/(sizeof(double)*(4*n_unique_alpha_+2*n_unique_beta_)))));
  if (!tangent && state.product_images) {
    state.product_images->extent=extent;
    state.product_images->states.resize(states);
    if (states==1) state.product_images->density=Eigen::MatrixXd::Zero(n,n);
  }
  const int panels=(n_unique_beta_+extent-1)/extent;
  const auto* support=product_support_.get();
  std::vector<int> reachable;
  std::vector<bool> reachable_mask(panels,false);
  if (support && !support->full()) {
    auto cursor=support->cursor(0); int b;
    while (cursor.next_column(b)) reachable_mask[b/extent]=true;
  } else std::fill(reachable_mask.begin(),reachable_mask.end(),true);
  for (int p=0;p<panels;++p) if (reachable_mask[p]) reachable.push_back(p);
  std::vector<std::weak_ptr<const detail::ColumnPattern>> patterns(panels);
  const auto pattern_at=[&](int panel) -> std::shared_ptr<const detail::ColumnPattern> {
    if (auto p=patterns[panel].lock()) return p;
    auto p=std::make_shared<detail::ColumnPattern>();
    const int count=std::min(extent,n_unique_beta_-panel*extent);
    p->full=!support || support->full(); p->offsets.resize(count+1); p->offsets[0]=0;
    std::vector<int> rows;
    for (int b=0;b<count;++b) {
      if (!p->full) {
        auto cursor=support->column(panel*extent+b,0); int a,beta; bool label;
        const auto begin=rows.size();
        while (cursor.next(a,beta,label)) rows.push_back(a);
        std::sort(rows.begin()+begin,rows.end());
      }
      p->offsets[b+1]=p->full ? (b+1)*n_unique_alpha_:rows.size();
    }
    p->rows=Eigen::Map<const Eigen::VectorXi>(rows.data(),rows.size());
    patterns[panel]=p; return p;
  };
  ExteriorTransformDirection da,db;
  if (delta_upper) { da=alpha.direction(*delta_upper); db=beta.direction(*delta_upper); }
  // Differentiated beta rows are shared by generation and left projection.
  detail::ColumnPanels beta_directions(cache_bytes,[&](int panel) {
    const int first=panel*extent,count=std::min(extent,n_unique_beta_-first);
    Eigen::MatrixXd rows=Eigen::MatrixXd::Zero(count,n_unique_beta_),dr=rows;
    rows.middleCols(first,count).setIdentity();
    beta.apply_directional_adjoint_right(db,&rows,&dr);
    auto p=std::make_shared<detail::ColumnPattern>();
    p->full=true; p->offsets.resize(count+1);
    for (int b=0;b<=count;++b) p->offsets[b]=b*n_unique_beta_;
    return detail::pack_columns(dr.transpose(),p);
  });
  const bool capture=diagnostic && diagnostic->capture_connection;
  if (capture) diagnostic->connection=StructureConnectionDirection{
      Eigen::MatrixXd::Zero(n_structures_,states),Eigen::MatrixXd::Zero(n_structures_,states),
      Eigen::MatrixXd::Zero(n_structures_,states),Eigen::VectorXd::Zero(states)};

  for (int state_id=0;state_id<states;++state_id) {
    const double weight=state.weights[state_id],energy=state.energies[state_id];
    if (weight==0 && !collect_images) continue;
    const bool cached=tangent && state.product_images && state.product_images->extent==extent &&
        !state.product_images->states[state_id].empty();
    const int sigma_width=cached ? 1:width;
    const std::size_t product_count=static_cast<std::size_t>(n_unique_alpha_)*n_unique_beta_;
    const bool resident_input=spin_products_.size()==product_count &&
        product_count<=cache_bytes/(width*sizeof(double));
    Eigen::MatrixXd primal_image,direction_image;
    std::size_t input_bytes=0;
    if (resident_input) {
      // Apply the ordered exterior jet once, instead of multiplying each beta
      // row panel by all original products. No orbital inverse is introduced.
      if (!cached || delta_upper)
        primal_image=expand_structure_block(state.structure_coefficients.col(state_id));
      if (tangent) direction_image=response ? expand_structure_block(response->col(state_id)):
          Eigen::MatrixXd::Zero(n_unique_alpha_,n_unique_beta_);
      if (delta_upper) {
        alpha.apply_directional_left(da,&primal_image,&direction_image);
        beta.apply_directional_right(db,&primal_image,&direction_image);
      } else {
        if (!cached) { alpha.apply_left(&primal_image); beta.apply_right(&primal_image); }
        if (tangent) { alpha.apply_left(&direction_image); beta.apply_right(&direction_image); }
      }
      input_bytes=sizeof(double)*(primal_image.size()+direction_image.size());
      if (cached) primal_image.resize(0,0);
    }
    Eigen::VectorXd coefficients;
    if (!resident_input && (!cached || delta_upper))
      coefficients=expand_support(state.structure_coefficients.col(state_id));
    Eigen::VectorXd q;
    if (response && !resident_input) q=expand_support(response->col(state_id));
    detail::ColumnPanels input(cache_bytes,[&](int panel) {
      ++result.images.tiles.input_panels;
      const int first=panel*extent,count=std::min(extent,n_unique_beta_-first);
      const auto rows=beta.row_block(first,count);
      Eigen::MatrixXd values=Eigen::MatrixXd::Zero(n_unique_alpha_,(cached ? 1:width)*count);
      Eigen::MatrixXd primal=Eigen::MatrixXd::Zero(n_unique_alpha_,count),dy=primal;
      Eigen::MatrixXd dr;
      if (delta_upper) dr=beta_directions.get(panel).unpack(n_unique_beta_).transpose();
      const int workers=coarse_openmp_thread_count(cooperative_openmp_thread_count(),
          static_cast<std::size_t>(n_unique_alpha_)*count);
      cooperative_openmp_for(n_unique_alpha_,workers,[&](std::ptrdiff_t a) {
        const auto first_product=mcscf::product_index(a,0,n_unique_beta_);
        auto group=std::lower_bound(spin_products_.begin(),spin_products_.end(),first_product);
        for (;group!=spin_products_.end() && *group<first_product+n_unique_beta_;++group) {
          const int id=group-spin_products_.begin(),b=*group-first_product;
          if (coefficients.size()) primal.row(a).noalias()+=coefficients[id]*rows->col(b).transpose();
          if (delta_upper) dy.row(a).noalias()+=coefficients[id]*dr.col(b).transpose();
          if (response) dy.row(a).noalias()+=q[id]*rows->col(b).transpose();
        }
      });
      if (delta_upper) alpha.apply_directional_left(da,&primal,&dy);
      else {
        if (!cached) {
          if (support) alpha.apply_supported_left(&primal,*support,first,count);
          else alpha.apply_left(&primal);
        }
        if (tangent) {
          if (support) alpha.apply_supported_left(&dy,*support,first,count);
          else alpha.apply_left(&dy);
        }
      }
      if (cached) values=std::move(dy);
      else {
        values.leftCols(count)=primal;
        if (tangent) values.rightCols(count)=dy;
      }
      return detail::pack_columns(std::move(values),pattern_at(panel));
    });
    // Independent readers borrow the same bounded primal/tangent panels.
    const auto reader=[&](int family) -> mcscf::DirectCiTileReader {
      return [&,family,position=std::vector<int>(n_unique_alpha_,-1)](
          const std::vector<int>& betas,const std::vector<int>& rows,
          Eigen::Ref<Eigen::MatrixXd> values) mutable {
        if (resident_input && (!cached || family==1)) {
          const auto& image=family==1 ? direction_image:primal_image;
          const int workers=coarse_openmp_thread_count(cooperative_openmp_thread_count(),
              rows.size()*betas.size());
          cooperative_openmp_for(betas.size(),workers,[&](std::ptrdiff_t col) {
            for (int a=0;a<static_cast<int>(rows.size());++a)
              values(a,col)=image(rows[a],betas[col]);
          });
          return;
        }
        for (int a=0;a<static_cast<int>(rows.size());++a) position[rows[a]]=a;
        const int component=cached ? (family==2 ? 1:0):family;
        const auto copy_column=[&](int col,const detail::ColumnTile& tile) {
          const int b=betas[col]%extent;
          const int first=tile.pattern->offsets[b],last=tile.pattern->offsets[b+1];
          if (tile.pattern->full) {
            for (int a=0;a<static_cast<int>(rows.size());++a)
              values(a,col)=tile.values(first+rows[a],component);
          } else {
            values.col(col).setZero();
            for (int k=first;k<last;++k) {
              const int a=position[tile.pattern->rows[k]];
              if (a>=0) values(a,col)=tile.values(k,component);
            }
          }
        };
        if (cached && family!=1) {
          // Immutable accepted panels span the entire requested column batch:
          // use one shared team, not one task launch per source panel.
          const int workers=coarse_openmp_thread_count(cooperative_openmp_thread_count(),
              rows.size()*betas.size());
          cooperative_openmp_for(betas.size(),workers,[&](std::ptrdiff_t col) {
            const int panel=betas[col]/extent;
            if (!reachable_mask[panel]) values.col(col).setZero();
            else copy_column(col,state.product_images->states[state_id][panel]);
          });
          for (const int a:rows) position[a]=-1;
          return;
        }
        for (int j=0;j<static_cast<int>(betas.size());) {
          const int panel=betas[j]/extent;
          int end=j+1;
          while (end<static_cast<int>(betas.size()) && betas[end]/extent==panel) ++end;
          if (!reachable_mask[panel]) { values.middleCols(j,end-j).setZero(); j=end; continue; }
          const auto& tile=input.get(panel);
          const int workers=bounded_openmp_thread_count(cooperative_openmp_thread_count(),end-j);
          cooperative_openmp_for(end-j,workers,[&](std::ptrdiff_t index) {
            copy_column(j+index,tile);
          });
          j=end;
        }
        for (const int a:rows) position[a]=-1;
      };
    };
    const auto y=reader(0),dy=reader(tangent ? 1:0);
    const auto column=[&](int family) -> mcscf::DirectCiColumnReader {
      return [&,family](int b) -> mcscf::DirectCiColumnView {
        if (!reachable_mask[b/extent]) return {nullptr,nullptr,0};
        if (resident_input && (!cached || family==1)) {
          const auto& image=family==1 ? direction_image:primal_image;
          return {image.col(b).data(),nullptr,n_unique_alpha_};
        }
        const auto& tile=cached && family!=1 ? state.product_images->states[state_id][b/extent]:input.get(b/extent);
        const int offset=tile.pattern->offsets[b%extent],last=tile.pattern->offsets[b%extent+1];
        return {tile.values.col(cached ? (family==2 ? 1:0):family).data()+offset,
            tile.pattern->full ? nullptr:tile.pattern->rows.data()+offset,last-offset};
      };
    };
    const auto y_column=column(0),dy_column=column(tangent ? 1:0);
    const auto residual_reader=reader(2);
    const auto residual_column=column(2);
    const mcscf::DirectCiTileReader read=[&](const std::vector<int>& betas,
        const std::vector<int>& rows,Eigen::Ref<Eigen::MatrixXd> values) {
      if (cached) { dy(betas,rows,values); return; }
      using View=Eigen::Map<Eigen::MatrixXd,0,Eigen::OuterStride<>>;
      View primal(values.data(),rows.size(),betas.size(),Eigen::OuterStride<>(values.outerStride()*width));
      y(betas,rows,primal);
      if (tangent) {
        View direction(values.data()+values.outerStride(),rows.size(),betas.size(),
            Eigen::OuterStride<>(values.outerStride()*width));
        dy(betas,rows,direction);
      }
    };
    const mcscf::DirectCiTileReader right=[&](const std::vector<int>& betas,
        const std::vector<int>& rows,Eigen::Ref<Eigen::MatrixXd> values) {
      values.setZero();
      const int offset=tangent && !cached ? 1:0;
      Eigen::Map<Eigen::MatrixXd,0,Eigen::OuterStride<>> view(
          values.data()+offset*values.outerStride(),rows.size(),betas.size(),
          Eigen::OuterStride<>(values.outerStride()*sigma_width));
      y(betas,rows,view);
    };
    // C^T is linear: first sum emitted panels in bounded source-group tiles,
    // then visit the structure incidence once per completed/evicted tile.
    // Neither the full touched support nor per-thread structure arrays are
    // required. Diagnostic outputs share the same cache and payload budget.
    const int projection_count=2+(capture ? 3:0);
    const std::size_t products=static_cast<std::size_t>(n_unique_alpha_)*n_unique_beta_;
    // Associativity removes the per-panel beta-row reconstruction on broad
    // support. Primal and left-direction images share the existing budget.
    const bool resident_projection=(collect_images || capture) &&
        spin_products_.size()==products && products<=cache_bytes/
            (sizeof(double)*projection_count*(delta_upper ? 2:1));
    std::vector<Eigen::MatrixXd> projected(projection_count),left_projected(projection_count);
    const auto group_extent=std::max<std::size_t>(1,std::min<std::size_t>(
        static_cast<std::size_t>(n_unique_alpha_)*extent,
        cache_bytes/(sizeof(double)*projection_count)));
    std::vector<double*> projection_targets;
    detail::ColumnPanels projections(cache_bytes,[&](int block) {
      const auto begin=static_cast<std::size_t>(block)*group_extent;
      const auto count=std::min(group_extent,spin_products_.size()-begin);
      auto pattern=std::make_shared<detail::ColumnPattern>();
      pattern->full=true; pattern->offsets.resize(2); pattern->offsets<<0,count;
      return detail::ColumnTile{pattern,Eigen::MatrixXd::Zero(count,projection_count)};
    },[&](int block,detail::ColumnTile& tile) {
      const auto begin=static_cast<std::size_t>(block)*group_extent;
      for (int target=0;target<static_cast<int>(projection_targets.size());++target)
        for (int local=0;local<tile.values.rows();++local) {
          const double value=tile.values(local,target);
          if (value==0.0) continue;
          const auto g=begin+local;
          for (auto t=spin_term_offsets_[g];t<spin_term_offsets_[g+1];++t) {
            const auto& term=spin_terms_[t];
            projection_targets[target][term.structure]+=term.coefficient*value;
          }
        }
    });
    const auto project=[&](int panel,const Eigen::MatrixXd& image,
        Eigen::Ref<Eigen::VectorXd> target,bool differentiate=false) {
      profiling::ScopedPhase phase(profiling::Phase::CiProductProjection);
      const int first=panel*extent,count=image.cols();
      Eigen::MatrixXd values=image,dvalues;
      if (differentiate) {
        dvalues=Eigen::MatrixXd::Zero(image.rows(),image.cols());
        alpha.apply_directional_adjoint_left(da,&values,&dvalues);
      } else if (support) alpha.apply_supported_left(&values,*support,first,count,true);
      else alpha.apply_adjoint_left(&values);
      auto found=std::find(projection_targets.begin(),projection_targets.end(),target.data());
      const int target_id=found-projection_targets.begin();
      if (found==projection_targets.end()) projection_targets.push_back(target.data());
      if (resident_projection) {
        if (!projected[target_id].size()) {
          projected[target_id]=Eigen::MatrixXd::Zero(n_unique_alpha_,n_unique_beta_);
          if (delta_upper) left_projected[target_id]=projected[target_id];
        }
        projected[target_id].middleCols(first,count)+=differentiate ? dvalues:values;
        if (differentiate) left_projected[target_id].middleCols(first,count)+=values;
        return;
      }
      const auto rows=beta.row_block(first,count);
      Eigen::MatrixXd dr;
      if (differentiate) dr=beta_directions.get(panel).unpack(n_unique_beta_).transpose();
      for (std::size_t begin=0;begin<spin_products_.size();begin+=group_extent) {
        auto& tile=projections.get(begin/group_extent);
        const auto size=std::min(group_extent,spin_products_.size()-begin);
        const int workers=coarse_openmp_thread_count(cooperative_openmp_thread_count(),size*count);
        cooperative_openmp_for(size,workers,[&](std::ptrdiff_t local) {
          const auto g=begin+local;
          const int a=spin_products_[g]/n_unique_beta_,b=spin_products_[g]%n_unique_beta_;
          tile.values(local,target_id)+=differentiate ? dvalues.row(a).dot(rows->col(b))+values.row(a).dot(dr.col(b)):
              values.row(a).dot(rows->col(b));
        });
      }
      ++result.images.tiles.structure_projections;
    };
    const auto sweep=[&]() {
      const int w=sigma_width;
      const bool forcing=capture && delta_sigma;
      const int output_width=w+(forcing ? 1:0);
      bool retain=false;
      detail::ColumnPanels output(cache_bytes,[&](int panel) {
        auto pattern=pattern_at(panel);
        return detail::ColumnTile{pattern,Eigen::MatrixXd::Zero(pattern->offsets[pattern->columns()],output_width)};
      },[&](int panel,detail::ColumnTile& tile) {
        ++result.images.tiles.hamiltonian_projections;
        const int first=panel*extent,count=tile.pattern->columns();
        Eigen::MatrixXd images=tile.unpack(n_unique_alpha_);
        const Eigen::MatrixXd h0=images.leftCols(count);
        if (weight!=0) {
          if (!tangent || cached)
            result.generator.noalias()+=2*weight*sigma.generator_adjoint_streamed(first,h0,
                {{y,y_column,1.0,false,cached}});
          else {
            result.generator.noalias()+=2*weight*sigma.generator_adjoint_streamed(first,h0,{{dy,dy_column}});
            result.generator.noalias()+=2*weight*sigma.generator_adjoint_streamed(
                first,images.middleCols((w-1)*count,count),{{y,y_column}});
          }
        }
        if (collect_images) {
          if (!tangent || cached) project(panel,h0,result.images.hamiltonian.col(state_id));
          else {
            project(panel,images.middleCols((w-1)*count,count),result.images.hamiltonian.col(state_id));
            if (delta_upper) project(panel,h0,result.images.hamiltonian.col(state_id),true);
          }
        }
        if (capture) {
          if (forcing) project(panel,images.middleCols(w*count,count),
              diagnostic->connection->integral_forcing.col(state_id));
          if (tangent) {
            project(panel,cached ? h0:Eigen::MatrixXd(images.middleCols((w-1)*count,count)),
                diagnostic->connection->representation_forcing.col(state_id));
            if (delta_upper && !cached)
              project(panel,h0,diagnostic->connection->representation_forcing.col(state_id),true);
          }
        }
        if (retain) {
          // Store Y and (H-E)Y on the exact support, sharing one pattern.
          // The next coefficient response then needs neither a new Y
          // transform nor a separate -E*D(Y,q) generator contraction.
          detail::ColumnTile accepted{tile.pattern,Eigen::MatrixXd(tile.values.rows(),2)};
          if (resident_input)
            accepted.values.col(0)=Eigen::Map<const Eigen::VectorXd>(
                primal_image.col(first).data(),accepted.values.rows());
          else accepted.values.col(0)=input.get(panel).values.col(0);
          accepted.values.col(1)=tile.values.col(0)-energy*accepted.values.col(0);
          state.product_images->bytes+=accepted.bytes();
          state.product_images->states[state_id][panel]=std::move(accepted);
        }
      });
      const auto writer=[&](int source_begin,int vectors,int destination_begin) -> mcscf::DirectCiTileWriter {
        return [&,source_begin,vectors,destination_begin](const std::vector<int>& betas,
          const std::vector<int>& rows,const mcscf::DirectCiColumnAccumulator& accumulate) {
        for (int j=0;j<static_cast<int>(betas.size());) {
          const int panel=betas[j]/extent;
          int end=j+1;
          while (end<static_cast<int>(betas.size()) && betas[end]/extent==panel) ++end;
          if (!reachable_mask[panel]) { j=end; continue; }
          auto& tile=output.get(panel);
          const int workers=cooperative_openmp_thread_count(),chunks=std::min<int>(workers,rows.size());
          const int columns=(end-j)*vectors;
          cooperative_openmp_for(columns*chunks,workers,[&](std::ptrdiff_t index) {
            const int column=index%columns,chunk=index/columns;
            const int beta=j+column/vectors,v=column%vectors,b=betas[beta]%extent;
            const int offset=tile.pattern->offsets[b],length=tile.pattern->offsets[b+1]-offset;
            if (length==0) return;
            accumulate(beta,source_begin+v,rows.size()*chunk/chunks,rows.size()*(chunk+1)/chunks,
                tile.values.col(destination_begin+v).data()+offset,
                tile.pattern->full ? nullptr:tile.pattern->rows.data()+offset,tile.pattern->full ? 0:length);
          });
          j=end;
        }
        };
      };
      const auto add=writer(0,w,0);
      {
        // Accepted bra and ket are identical. Keep the empty-reader identity
        // so the channel kernel reuses its scatter; apply the weight outside.
        const auto gradient=delta_sigma ? sigma.apply_tiled_direction(w,read,right,*delta_sigma,add,
            forcing ? writer(w-1,1,w):mcscf::DirectCiTileWriter{},support,spin_exchange_parity_):
            sigma.apply_tiled_with_adjoint(w,read,add,
                tangent ? right:mcscf::DirectCiTileReader{},support,spin_exchange_parity_);
        // Keep state weights outside the shared scatter and both kernel actions.
        result.one_electron.noalias()+=(tangent ? 2*weight:weight)*gradient.one_electron;
        result.pair_kernel.noalias()+=(tangent ? 2*weight:weight)*gradient.pair_kernel;
        ++result.images.tiles.integral_adjoint_sweeps;
      }
      if (!tangent && state.product_images && output.evictions==0)
        for (const int panel:reachable) output.get(panel);
      if (!tangent && state.product_images && output.evictions==0 &&
          2*output.peak_bytes<=cache_bytes-state.product_images->bytes) {
        retain=true;
        state.product_images->states[state_id].resize(panels);
      }
      output.finish();
      result.images.tiles.output_evictions+=output.evictions;
      result.images.tiles.output_resident_bytes=std::max(result.images.tiles.output_resident_bytes,output.peak_bytes);
    };
    if (cached && collect_images && delta_upper) {
      // Reuse accepted Y/residual panels for all left-direction/backward consumers.
      // No per-HVP replay of the accepted sigma is needed when it is resident.
      const auto& images=state.product_images->states[state_id];
      for (const int panel:reachable) {
        if (!images[panel].pattern) continue;
        const Eigen::MatrixXd values=images[panel].unpack(n_unique_alpha_);
        const int count=images[panel].pattern->columns();
        const Eigen::MatrixXd residual=values.rightCols(count);
        {
          const Eigen::MatrixXd h=residual+energy*values.leftCols(count);
          project(panel,h,result.images.hamiltonian.col(state_id),true);
          if (capture) project(panel,h,diagnostic->connection->representation_forcing.col(state_id),true);
        }
      }
    }
    sweep();
    double norm=0,parallel=0,square=0;
    std::vector<int> all_rows(n_unique_alpha_); std::iota(all_rows.begin(),all_rows.end(),0);
    for (const int panel:reachable) {
      const int first=panel*extent,count=std::min(extent,n_unique_beta_-first);
      std::vector<int> betas(count); std::iota(betas.begin(),betas.end(),first);
      Eigen::MatrixXd values,dv;
      if (!tangent || delta_upper || capture) {
        values.resize(n_unique_alpha_,count);
        y(betas,all_rows,values);
      }
      if (tangent) { dv.resize(n_unique_alpha_,count); dy(betas,all_rows,dv); }
      if (weight!=0) {
        if (tangent) {
          // D(R,Z)=D(Z,R)^T: both accepted cotangents share the generated
          // Z panel and one spin-link traversal, without unpacking R first.
          if (cached) result.generator.noalias()+=2*weight*sigma.generator_adjoint_streamed(
              first,dv,{{y,y_column,-energy,false,true},{residual_reader,residual_column,1.0,true,true}});
          else result.generator.noalias()-=2*weight*energy*
              sigma.generator_adjoint_streamed(first,dv,{{y,y_column}},true);
        }
        else {
          const auto density=sigma.generator_adjoint_streamed(first,values,{{y,y_column}});
          result.generator.noalias()-=2*weight*energy*density;
          if (states==1 && state.product_images) state.product_images->density+=density;
        }
      }
      if (collect_images) {
        if (!tangent) project(panel,values,result.images.overlap.col(state_id));
        else {
          project(panel,dv,result.images.overlap.col(state_id));
          if (delta_upper) project(panel,values,result.images.overlap.col(state_id),true);
        }
      }
      if (capture) {
        project(panel,dv,diagnostic->connection->overlap.col(state_id));
        norm+=values.squaredNorm(); parallel+=(values.array()*dv.array()).sum(); square+=dv.squaredNorm();
      }
    }
    projections.finish();
    if (resident_projection) {
      profiling::ScopedPhase phase(profiling::Phase::CiProductProjection);
      for (int target=0;target<static_cast<int>(projection_targets.size());++target) {
        if (delta_upper)
          beta.apply_directional_adjoint_right(db,&left_projected[target],&projected[target]);
        else beta.apply_adjoint_right(&projected[target]);
        Eigen::Map<Eigen::VectorXd>(projection_targets[target],n_structures_)+=
            contract_spin_product_block(projected[target]).col(0);
        ++result.images.tiles.structure_projections;
      }
    }
    if (capture) {
      if (delta_sigma) diagnostic->connection->representation_forcing.col(state_id)-=
          diagnostic->connection->integral_forcing.col(state_id);
      diagnostic->connection->squared_norm[state_id]=square-(norm>0 ? parallel*parallel/norm:0);
      diagnostic->connection->representation_forcing.col(state_id).noalias()-=energy*result.images.overlap.col(state_id);
    }
    result.images.tiles.input_evictions+=input.evictions;
    result.images.tiles.input_resident_bytes=std::max(result.images.tiles.input_resident_bytes,input.peak_bytes);
    result.images.tiles.input_resident_bytes=std::max(result.images.tiles.input_resident_bytes,input_bytes);
  }
  if (!tangent && states==1 && state.product_images)
    state.product_images->density=state.product_images->density.selfadjointView<Eigen::Upper>();
  return result;
}

} // namespace xmvb::vb
