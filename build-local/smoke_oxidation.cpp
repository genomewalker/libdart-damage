#include <iostream>
#include <cmath>
#include "taph/library_interpretation.hpp"
#include "taph/sample_damage_profile.hpp"

int main(){
  taph::SampleDamageProfile dp;
  dp.n_reads = 200000;
  dp.d_max_5prime = 0.12f;
  dp.d_max_3prime = 0.04f;
  dp.terminal_z_5prime = 9.0f;
  dp.terminal_z_3prime = 5.0f;
  dp.mixture_converged = true;
  dp.mixture_identifiable = true;
  dp.mixture_d_ancient = 0.11f;
  dp.mixture_pi_ancient = 0.42f;

  dp.baseline_g_to_t_count = 1200;
  dp.baseline_g_total = 40000;
  dp.baseline_c_to_a_count = 900;
  dp.baseline_c_ox_total = 35000;

  dp.baseline_a_to_t_count = 300;
  dp.baseline_a_total = 42000;
  dp.baseline_c_to_g_count = 260;
  dp.baseline_c_bg_total = 36000;

  dp.ox_gt_rate_terminal = 0.041f;
  dp.ox_ca_rate_terminal = 0.038f;
  dp.ox_gt_rate_interior = 0.028f;
  dp.ox_ca_rate_interior = 0.026f;

  auto r = taph::compute_preservation_summary(dp,false,false,false,2.0,1.5,0.8,0.5);

  if(!std::isfinite(r.oxidation.raw_rate.estimate)) return 2;
  if(!std::isfinite(r.oxidation.bg_rate.estimate)) return 3;
  if(!std::isfinite(r.oxidation.excess_rate.estimate)) return 4;
  if(!std::isfinite(r.oxidation.reliability_score)) return 5;
  std::cout << "raw=" << r.oxidation.raw_rate.estimate
            << " bg=" << r.oxidation.bg_rate.estimate
            << " excess=" << r.oxidation.excess_rate.estimate
            << " rel=" << r.oxidation.reliability
            << " rel_score=" << r.oxidation.reliability_score
            << " term_enrich=" << r.oxidation.terminal_enrichment
            << " legacy_eff=" << r.oxidation_eff
            << " legacy_ev=" << r.oxidation_evidence
            << "\n";
  return 0;
}
