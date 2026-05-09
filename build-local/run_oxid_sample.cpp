#include <iostream>
#include <string>
#include "taph/frame_selector_decl.hpp"
#include "taph/library_interpretation.hpp"

int main() {
    taph::SampleDamageProfile dp;
    std::string seq;
    size_t n = 0;
    while (std::getline(std::cin, seq)) {
        if (seq.empty()) continue;
        taph::FrameSelector::update_sample_profile(dp, seq);
        ++n;
    }
    taph::FrameSelector::finalize_sample_profile(dp);
    const bool is_ss = (dp.library_type == taph::SampleDamageProfile::LibraryType::SINGLE_STRANDED);
    auto cpg = taph::compute_cpg_score(dp);
    auto oxi = taph::compute_oxog_interior_score(dp);
    auto otr = taph::compute_oxog_trinuc(dp);
    auto pres = taph::compute_preservation_summary(dp, is_ss, false, false, cpg.z, oxi.z, otr.cosine, 1.0);

    std::cout << "reads=" << n << " lib=" << dp.library_type_str() << "\n";
    std::cout << "raw=" << pres.oxidation.raw_rate.estimate
              << " raw_ci=[" << pres.oxidation.raw_rate.ci95_low << "," << pres.oxidation.raw_rate.ci95_high << "]\n";
    std::cout << "control=" << pres.oxidation.control_rate.estimate
              << " control_ci=[" << pres.oxidation.control_rate.ci95_low << "," << pres.oxidation.control_rate.ci95_high << "]\n";
    std::cout << "excess=" << pres.oxidation.excess_rate.estimate
              << " excess_ci=[" << pres.oxidation.excess_rate.ci95_low << "," << pres.oxidation.excess_rate.ci95_high << "]\n";
    std::cout << "z=" << pres.oxidation.z_score
              << " bins=" << pres.oxidation.bins_used
              << " eff_bins=" << pres.oxidation.effective_bins
              << " heterogeneity=" << pres.oxidation.heterogeneity
              << " reliability=" << pres.oxidation.reliability
              << " rel_score=" << pres.oxidation.reliability_score << "\n";
    std::cout << "legacy_eff=" << pres.oxidation_eff
              << " legacy_evidence=" << pres.oxidation_evidence << "\n";
    return 0;
}
