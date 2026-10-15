/** @file Exact-support sparse exterior transform and sampled CI connectivity. */
#include <algorithm>
#include <bit>
#include <chrono>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <random>
#include <stdexcept>
#include <vector>

using Clock = std::chrono::steady_clock;
struct Entry { std::uint64_t key; double value; };
static_assert(sizeof(Entry)==16);
double seconds(Clock::time_point t) {
  return std::chrono::duration<double>(Clock::now()-t).count();
}

/** Apply one spin shear by merging two ordered supports; never screen values. */
void shear(std::vector<Entry>& x, int p, int q, double t,
           std::size_t& peak_workspace) {
  if (t==0.) return;
  const auto bp=std::uint64_t{1}<<p, bq=std::uint64_t{1}<<q;
  const auto between=((bq-1)^(bp-1))^bp;
  std::vector<Entry> added;
  for (const auto& e:x) {
    if (!(e.key&bq) || (e.key&bp)) continue;
    const double phase=(std::popcount(e.key&between)&1)?-1.:1.;
    const double value=t*phase*e.value;
    if (value!=0.) added.push_back({e.key^bp^bq,value});
  }
  if (added.empty()) return;
  // Replacement shifts every eligible key by the same integer: order persists.
  std::vector<Entry> result;
  result.reserve(x.size()+added.size());
  std::size_t i=0,j=0;
  while(i<x.size() || j<added.size()) {
    if(j==added.size() || (i<x.size() && x[i].key<added[j].key))
      result.push_back(x[i++]);
    else if(i==x.size() || added[j].key<x[i].key)
      result.push_back(added[j++]);
    else {
      const auto key=x[i].key;
      const double value=x[i++].value+added[j++].value;
      if(value!=0.) result.push_back({key,value});
    }
  }
  peak_workspace=std::max(peak_workspace,
      sizeof(Entry)*(x.capacity()+added.capacity()+result.capacity()));
  x.swap(result);
}

/** Generate all spin-preserving single/double Slater-Condon candidates. */
std::uint64_t count_hits(const std::vector<Entry>& x, std::uint64_t key,int n) {
  std::vector<std::uint64_t> singles[2], doubles[2];
  for(int spin=0;spin<2;++spin) {
    std::vector<int> occ,virt;
    for(int p=spin*n;p<(spin+1)*n;++p)
      ((key>>p)&1?occ:virt).push_back(p);
    for(int i:occ) for(int a:virt)
      singles[spin].push_back((1ULL<<i)^(1ULL<<a));
    for(std::size_t i=0;i<occ.size();++i)
      for(std::size_t j=0;j<i;++j)
        for(std::size_t a=0;a<virt.size();++a)
          for(std::size_t b=0;b<a;++b)
            doubles[spin].push_back((1ULL<<occ[i])^(1ULL<<occ[j])^
                                    (1ULL<<virt[a])^(1ULL<<virt[b]));
  }
  std::uint64_t hits=1;
  const auto visit=[&](std::uint64_t delta) {
    const auto target=key^delta;
    auto it=std::lower_bound(x.begin(),x.end(),target,
        [](const Entry& e,std::uint64_t k){return e.key<k;});
    hits+=(it!=x.end() && it->key==target);
  };
  for(int spin=0;spin<2;++spin) {
    for(auto delta:singles[spin]) visit(delta);
    for(auto delta:doubles[spin]) visit(delta);
  }
  for(auto a:singles[0]) for(auto b:singles[1]) visit(a^b);
  return hits;
}

int main() {
  const auto start=Clock::now();
  int n;
  std::ifstream input("orbitals.txt");
  input>>n;
  if(!input || n>30 || n%2) throw std::runtime_error("invalid orbital input");
  std::vector<double> r(n*n);
  for(auto& value:r) input>>value;
  if(!input) throw std::runtime_error("incomplete orbital input");
  const int m=n/2;
  std::vector<Entry> x;
  for(std::uint64_t choice=0;choice<(1ULL<<m);++choice) {
    std::uint64_t a=0;
    for(int g=0;g<m;++g) a|=1ULL<<(2*g+((choice>>g)&1));
    const auto b=((1ULL<<n)-1)^a;
    double value=std::pow(2.,-m/2.);
    for(int p=0;p<n;++p) value*=r[p*n+p];
    x.push_back({a|(b<<n),value});
  }
  std::sort(x.begin(),x.end(),[](const Entry&a,const Entry&b){return a.key<b.key;});
  std::vector<double> unit=r;
  for(int j=0;j<n;++j) for(int i=0;i<n;++i) unit[i*n+j]/=r[j*n+j];
  std::size_t peak=x.size(),peak_workspace=x.capacity()*sizeof(Entry);
  std::ofstream trace("support.csv");
  trace<<"removed,inserted,spin,support,seconds\n";
  for(int q=1;q<n;++q) for(int p=q-1;p>=0;--p) {
    const double t=unit[p*n+q];
    for(int row=0;row<n;++row) unit[row*n+q]-=t*unit[row*n+p];
    for(int spin=0;spin<2;++spin) {
      shear(x,p+spin*n,q+spin*n,t,peak_workspace);
      peak=std::max(peak,x.size());
      trace<<q<<','<<p<<','<<spin<<','<<x.size()<<','<<seconds(start)<<'\n';
    }
  }
  trace.close();
  const double transform_seconds=seconds(start);
  long double norm=0.;
  std::size_t significant=0;
  for(auto e:x) {norm+=static_cast<long double>(e.value)*e.value;
                significant+=std::abs(e.value)>1e-14;}
  double reference_error=0.;
  std::ifstream reference("reference.bin",std::ios::binary);
  if(reference) {
    Entry e;
    std::size_t i=0;
    while(reference.read(reinterpret_cast<char*>(&e),sizeof(e))) {
      while(i<x.size() && x[i].key<e.key) {
        reference_error=std::max(reference_error,std::abs(x[i++].value));
      }
      const double value=(i<x.size() && x[i].key==e.key)?x[i++].value:0.;
      reference_error=std::max(reference_error,std::abs(value-e.value));
    }
    while(i<x.size()) reference_error=std::max(reference_error,std::abs(x[i++].value));
    if(reference_error>1e-11) throw std::runtime_error("exterior reference mismatch");
  }
  std::ofstream result("result.json");
  result<<std::setprecision(17)<<"{\n\"n\":"<<n<<",\n\"support\":"<<x.size()
        <<",\n\"peak_support\":"<<peak<<",\n\"significant_above_1e_14\":"<<significant
        <<",\n\"norm\":"<<static_cast<double>(norm)
        <<",\n\"transform_seconds\":"<<transform_seconds
        <<",\n\"coefficient_bytes\":"<<x.size()*sizeof(Entry)
        <<",\n\"peak_vector_workspace_bytes\":"<<peak_workspace
        <<",\n\"reference_max_error\":"<<reference_error<<"\n}\n";
  result.close();
  std::cout<<"transform completed n="<<n<<" support="<<x.size()
           <<" seconds="<<transform_seconds<<std::endl;
  const auto count_start=Clock::now();
  const std::size_t samples=std::min<std::size_t>(4096,x.size());
  std::mt19937_64 random(20261015);
  std::uniform_int_distribution<std::size_t> pick(0,x.size()-1);
  long double sum=0.,sum2=0.;
  for(std::size_t row=0;row<samples;++row) {
    const auto index=samples==x.size()?row:pick(random);
    const auto hits=count_hits(x,x[index].key,n);
    sum+=hits;sum2+=static_cast<long double>(hits)*hits;
  }
  const double mean=static_cast<double>(sum/samples);
  const double se=samples==x.size()?0.:
      std::sqrt(std::max(0.,static_cast<double>((sum2-sum*sum/samples)/(samples-1)/samples)));
  std::ofstream counts("connectivity.json");
  counts<<std::setprecision(17)<<"{\n\"sampled_rows\":"<<samples
        <<",\n\"mean_potential_connections_in_support\":"<<mean
        <<",\n\"mean_standard_error\":"<<se
        <<",\n\"estimated_directed_connections\":"<<mean*x.size()
        <<",\n\"seconds\":"<<seconds(count_start)<<"\n}\n";
}
