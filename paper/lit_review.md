# Literature Review Working Notes

55 sources across 5 buckets. Citation keys match `refs.bib`.

## Buckets

- **A.** Migraine GWAS and genetics (12)
- **B.** CGRP biology and pharmacology (13)
- **C.** Autonomic / POTS biology and comorbidity (12)
- **D.** Computational pathway and network methods (12)
- **E.** Existing triangulation work — gap analysis (6 unique to E; A and D entries also load-bear)

---

## A. Migraine GWAS and genetics

### `hautakangas_2022_migraine_gwas` — anchor paper
Hautakangas et al. 2022, *Nat Genet* 54(2):152–160. doi:10.1038/s41588-021-00990-0.
**Key claim**: Largest published migraine GWAS (102,084 cases / 771,257 controls); 123 risk loci, 86 novel. Tissue/cell-type enrichment in vascular and CNS tissues. New chromosome 11 locus contains *CALCA*/*CALCB* (CGRP ligands), but the paper explicitly notes *CALCRL*, *RAMP1*, and *RCP* (receptor genes) "show no statistically comparable association (all P > 10⁻⁴)." Subtype-specific loci: *HMOX2*, *CACNA1A*, *MPPED2* (with aura); *SPINK2*, *FECH* (without aura). MAGMA used MSigDB curated + GO collections; **no autonomic/adrenergic-specific gene set tested**.
**Novelty implication**: The CGRP **ligand vs receptor asymmetry** is the load-bearing fact for our triangulation. Generic enrichment exists; targeted autonomic/adrenergic enrichment does not.

### `gormley_2016_migraine_gwas`
Gormley et al. 2016, *Nat Genet* 48(8):856–866. doi:10.1038/ng.3598.
**Key claim**: Pre-Hautakangas anchor; 38 loci, 28 novel (60K cases). Tissue/cell-type enrichment highlighted vascular and smooth-muscle expression. CGRP-pathway genes were not flagged.
**Novelty implication**: Sets historical baseline; no autonomic gene-set analysis.

### `anttila_2013_migraine_gwas`
Anttila et al. 2013, *Nat Genet* 45(8):912–917. doi:10.1038/ng.2676.
**Key claim**: Earliest large meta-GWAS (23K cases); 12 loci. No CGRP loci.
**Novelty implication**: Historical context; lapped 10× by Hautakangas.

### `choquet_2021_multiethnic`
Choquet et al. 2021, *Commun Biol* 4:864. doi:10.1038/s42003-021-02356-y.
**Key claim**: Multiethnic meta-analysis; 79 loci, 45 novel. Independently corroborates *CALCB* migraine locus. Female-specific loci: *CPS1*, *PBRM1*, *SLC25A21*. FUMA enrichment surfaced flavonoid glucuronide / ascorbate metabolism — **no autonomic gene-set test**.
**Novelty implication**: Strengthens CGRP-ligand evidence; sex-specific loci hint at hormonal axis.

### `bjornsdottir_2023_rare_migraine_variants`
Bjornsdottir et al. 2023, *Nat Genet* 55(11):1843–1853. doi:10.1038/s41588-023-01538-0.
**Key claim**: deCODE rare-variant analysis (~1.3M individuals). Frameshift in *PRRT2* confers large migraine-with-aura risk (OR≈5.4); *SCN11A* LoF protective; *KCNK5* cis-regulatory protective. Aura vs no-aura genetic divergence at rare-variant tier. **No CGRP-pathway rare variants reported.**
**Novelty implication**: Aura/no-aura divergence relevant if our triangulation distinguishes POTS-comorbid migraine subtypes.

### `hautakangas_2025_migraine_finemapping`
Hautakangas et al. 2025, *Nat Commun* 16. doi:10.1038/s41467-025-64880-3.
**Key claim**: Fine-mapping of 98K cases; 181 credible sets; 7 variants with PIP > 0.9. CGRP-pathway genes are **not** credible-set drivers in main results.
**Novelty implication**: Reinforces ligand/receptor asymmetry; no autonomic gene-set analysis.

### `ghaffar_2023_migraine_eqtl_twas`
Ghaffar & Nyholt 2023, *Hum Genet* 142(8):1113–1137. doi:10.1007/s00439-023-02568-8.
**Key claim**: TWAS over GTEx v8 / 49 tissues; 62 novel candidate genes. CGRP-pathway genes not in TWAS hits; tissue prioritization vascular/inflammatory; **no autonomic enrichment**.
**Novelty implication**: TWAS-style integration has not surfaced CGRP receptor or autonomic genes.

### `liu_2024_han_chinese_migraine_gwas`
Liu et al. 2024, *J Clin Neurol* 20(4):439–449. doi:10.3988/jcn.2023.0331.
**Key claim**: Han Chinese GWAS (1,561 cases). Novel loci *DDX1/LINC01804*, *ELMO1*. Underpowered relative to European mega-analyses.
**Novelty implication**: Trans-ancestry generalizability; no autonomic gene-set test.

### `alfayyadh_2024_hemiplegic_migraine_review`
Alfayyadh et al. 2024, *Genes* 15(4):443. doi:10.3390/genes15040443.
**Key claim**: Reviews monogenic hemiplegic migraine: *CACNA1A* (FHM1), *ATP1A2* (FHM2), *SCN1A* (FHM3) on glutamate-clearance / cortical-spreading-depression susceptibility. ~75% of HM patients lack mutations in these three. **Does not address autonomic mechanisms.**
**Novelty implication**: Anchors monogenic FHM context; explicit gap on autonomic dimension.

### `staehr_2025_fhm_ukb_exome`
Staehr et al. 2025, *Cephalalgia* 45(1). doi:10.1177/03331024241306103.
**Key claim**: UK Biobank exome-wide association: heterozygous FHM-gene variants confer sub-Mendelian migraine risk in population biobanks.
**Novelty implication**: The same population-biobank exome logic could apply to CGRP-pathway and autonomic genes — has not been done.

### `olofsson_2024_migraine_twin_review`
Olofsson 2024, *Headache* 64(8):1049–1058. doi:10.1111/head.14789.
**Key claim**: Adult migraine heritability concentrates at 0.36–0.48 across twin studies.
**Novelty implication**: Common-variant GWAS still leaves substantial missing heritability — motivates pathway-level integration.

### `sutherland_2024_migraine_genetics_review`
Sutherland et al. 2024, *Lancet Neurol* 23(4):429–446. doi:10.1016/S1474-4422(24)00026-7.
**Key claim**: 2024 narrative review framing migraine genetics as excitatory–inhibitory imbalance + neurovascular convergence. Discusses CGRP as vasoactive neuropeptide but **does not synthesize CGRP-pathway gene-set evidence**.
**Novelty implication**: State-of-field review's silence on autonomic gene-set analyses confirms our gap.

### Bucket A synthesis
Migraine GWAS is mature (Anttila 2013 → Gormley 2016 → Hautakangas 2022 → Hautakangas 2025 fine-mapping) with multiethnic confirmation (Choquet 2021) and rare-variant frontier opening (Bjornsdottir 2023). Twin heritability is ~0.36–0.48. **The CGRP ligand locus (*CALCA*/*CALCB*) is genome-wide significant; the receptor-complex genes (*CALCRL*, *RAMP1*, *RCP*) are not.** No major migraine GWAS has run pathway enrichment specifically on autonomic/adrenergic gene sets — Hautakangas 2022 used MSigDB+GO via MAGMA, Gormley 2016 reported broad vascular/smooth-muscle tissue enrichment, Choquet 2021's FUMA surfaced metabolic pathways, Ghaffar 2023's TWAS prioritized vascular/inflammatory tissues, and the 2024 Sutherland review does not catalogue any autonomic analysis.

---

## B. CGRP biology and pharmacology

### `amara_1982_calca`
Amara et al. 1982, *Nature* 298(5871):240–244. doi:10.1038/298240a0.
**Key claim**: Foundational paper establishing tissue-specific alternative RNA processing of CALCA — calcitonin mRNA in thyroid C-cells, CGRP-encoding mRNA in neurons.
**Novelty implication**: Anchors *CALCA* biology in framing.

### `liang_2018_cgrpr_cryoem`
Liang et al. 2018, *Nature* 561(7724):492–497. doi:10.1038/s41586-018-0535-y.
**Key claim**: 3.3 Å cryo-EM of active human CGRP receptor (CLR + RAMP1) with CGRP and Gs (PDB 6E3Y). RAMP1 sits at CLR TM3/4/5 interface; minimal direct CGRP contact, consistent with allosteric modulation.
**Novelty implication**: Provides molecular structure for any computational structural work.

### `lassen_2002_cgrp_provocation`
Lassen et al. 2002, *Cephalalgia* 22(1):54–61. doi:10.1046/j.1468-2982.2002.00310.x.
**Key claim**: Double-blind crossover (n=12): IV α-CGRP provoked headache in all 12 migraineurs without aura vs 1/12 placebo. First human evidence CGRP is sufficient to provoke migraine-like attacks.
**Novelty implication**: Central perturbational evidence for CGRP causality.

### `goadsby_2017_strive_erenumab`
Goadsby et al. 2017, *NEJM* 377(22):2123–2132. doi:10.1056/NEJMoa1705848.
**Key claim**: Phase 3 STRIVE (n=955): erenumab (anti-CGRP-receptor mAb) reduced monthly migraine days by 3.2 (70mg) / 3.7 (140mg) vs 1.8 placebo over 6 months. ≥50% responder rates 43%/50% vs 27%.
**Novelty implication**: Anti-CGRP-receptor effect-size benchmark.

### `silberstein_2017_halo_fremanezumab`
Silberstein et al. 2017, *NEJM* 377(22):2113–2122. doi:10.1056/NEJMoa1709038.
**Key claim**: HALO-CM (n=1130 chronic migraine): fremanezumab (anti-CGRP ligand mAb) reduced monthly headache days by 4.6/4.3 vs 2.5 placebo. ≥50% responder rates 41%/38% vs 18%.
**Novelty implication**: Chronic-migraine ligand-targeting mAb anchor.

### `stauffer_2018_evolve1_galcanezumab`
Stauffer et al. 2018, *JAMA Neurol* 75(9):1080–1088. doi:10.1001/jamaneurol.2018.1212.
**Key claim**: EVOLVE-1 (n=858): galcanezumab reduced monthly migraine days by ~4.7 vs ~2.8 placebo. ≥50% responder rates ~62% vs 39%.
**Novelty implication**: Episodic-migraine ligand-targeting effect size.

### `lipton_2020_promise2_eptinezumab`
Lipton et al. 2020, *Neurology* 94(13):e1365–e1377. doi:10.1212/WNL.0000000000009169.
**Key claim**: PROMISE-2 (n=1072): IV eptinezumab reduced monthly migraine days by ~7.7/8.2 vs 5.6 placebo. Day-1 post-infusion migraine probability dropped 42% → 28%.
**Novelty implication**: Rapid IV onset demonstrates pharmacological precision.

### `croop_2019_rimegepant`
Croop et al. 2019, *Lancet* 394(10200):737–745. doi:10.1016/S0140-6736(19)31606-X.
**Key claim**: Phase 3 acute migraine, oral rimegepant 75mg ODT: 2h pain freedom 21% vs 11% placebo. First gepant pivotal trial.
**Novelty implication**: Small-molecule CGRP receptor antagonism for acute treatment.

### `dodick_2019_achieve1_ubrogepant`
Dodick et al. 2019, *NEJM* 381(23):2230–2241. doi:10.1056/NEJMoa1813049.
**Key claim**: ACHIEVE-I (n=1672): ubrogepant 50mg/100mg achieved 2h pain freedom 19%/21% vs 12% placebo.
**Novelty implication**: Second pivotal acute gepant.

### `ailani_2021_advance_atogepant`
Ailani et al. 2021, *NEJM* 385(8):695–706. doi:10.1056/NEJMoa2035908.
**Key claim**: ADVANCE (n=910): daily oral atogepant 10/30/60mg reduced monthly migraine days by 3.7/3.9/4.2 vs 2.5 placebo. First daily oral gepant for prevention.
**Novelty implication**: Small-molecule CGRP antagonism is sufficient for prophylaxis.

### `edvinsson_2018_cgrp_review`
Edvinsson et al. 2018, *Nat Rev Neurol* 14(6):338–350. doi:10.1038/s41582-018-0003-1.
**Key claim**: Authoritative CGRP-in-migraine review; argues trigeminal ganglion (peripheral, outside BBB) is likely primary site of action.
**Novelty implication**: Primary review citation for CGRP-migraine mechanism.

### `russell_2014_cgrp_physiology`
Russell et al. 2014, *Physiol Rev* 94(4):1099–1142. doi:10.1152/physrev.00034.2013.
**Key claim**: Comprehensive CGRP physiology review: CGRP is among the most potent endogenous vasodilators, mediates neurogenic inflammation, and counteracts the renin–angiotensin–aldosterone and sympathetic systems in cardiovascular regulation.
**Novelty implication**: Directly supports our triangulation thesis — CGRP is a sensory–autonomic–vascular peptide, not solely a pain peptide.

### `saely_2021_erenumab_hypertension`
Saely et al. 2021, *Headache* 61(1):202–208. doi:10.1111/head.14051.
**Key claim**: FAERS postmarketing case series: 61 cases of elevated BP on erenumab (median +39/+28 mmHg). Led to 2020 US label update adding hypertension to Warnings.
**Novelty implication**: Caveat for Limitations; **also corroborates autonomic framing** — if blocking CGRP raises BP, endogenous CGRP is doing autonomic work at baseline.

### Bucket B synthesis
The molecular case for CGRP rests on three converging lines: mechanistic (Amara 1982 splicing; Liang 2018 cryo-EM), perturbational (Lassen 2002 infusion provokes attacks), and pharmacological (four pivotal mAb trials with consistent ~3.5–4.7 monthly migraine day reductions and ~40–60% responder rates; three pivotal gepant trials). Edvinsson 2018 consolidates these into the trigeminovascular framework. Critically, Russell 2014 documents CGRP's role as a tonic counter-regulator of sympathetic and renin–angiotensin–aldosterone tone, placing it squarely in autonomic-vascular territory. Saely 2021 is doubly informative: it constrains Limitations (anti-CGRP receptor blockade carries a real BP signal) and corroborates autonomic framing.

---

## C. Autonomic / POTS biology and comorbidity

### `sheldon_2015_hrs`
Sheldon et al. 2015, *Heart Rhythm* 12(6):e41–e63. doi:10.1016/j.hrthm.2015.03.029.
**Key claim**: HRS consensus: POTS = sustained HR increase ≥30 bpm (≥40 in adolescents) within 10 min of standing, no orthostatic hypotension, chronic symptoms ≥3 months.
**Novelty implication**: Anchors POTS phenotype definition.

### `vernino_2021_nih`
Vernino et al. 2021, *Auton Neurosci* 235:102828. doi:10.1016/j.autneu.2021.102828.
**Key claim**: NIH Expert Consensus: POTS is heterogeneous (neuropathic / hyperadrenergic / hypovolemic mechanisms overlap), 80–94% female, principal comorbidities are migraine, EDS, MCAS, ME/CFS. **POTS pathophysiology is multifactorial — no unified mechanism.**
**Novelty implication**: Most authoritative recent statement of the open mechanistic territory we enter.

### `raj_2013_pots`
Raj 2013, *Circulation* 127(23):2336–2342. doi:10.1161/CIRCULATIONAHA.112.144501.
**Key claim**: Subtype taxonomy — neuropathic (lower-limb sympathetic denervation), hyperadrenergic (upright NE >600 pg/mL with rising BP), hypovolemic (reduced RBC + plasma volume); subtypes overlap in most patients.
**Novelty implication**: Canonical subtype taxonomy.

### `shannon_2000_net`
Shannon et al. 2000, *NEJM* 342(8):541–549. doi:10.1056/NEJM200002243420803.
**Key claim**: First identified monogenic cause of orthostatic tachycardia: heterozygous *SLC6A2* (NET) A457P with >98% LoF in a family with familial orthostatic intolerance.
**Novelty implication**: Signature monogenic precedent — autonomic gene defect → POTS-like phenotype.

### `nakao_2012_gnb3`
Nakao et al. 2012, *Pediatr Int* 54(6):829–837. doi:10.1111/j.1442-200X.2012.03707.x.
**Key claim**: *GNB3* C825T TT genotype overrepresented in pediatric POTS (45.8% vs 20.0%, p=0.036).
**Novelty implication**: Rare candidate-gene replication signal in POTS.

### `qu_2025_pots_genetics` ⭐
Qu et al. 2025, *Clin Auton Res* 35(3):431–451. doi:10.1007/s10286-025-01110-2.
**Key claim**: First systematic GWAS+WES of POTS (207 European cases / 4063 controls; WES on 87). No common variants reached genome-wide significance. Over-representation analysis on the 5,670 nominally-significant SNPs (P<0.05) recovered HALLMARK_ESTROGEN_RESPONSE_EARLY (P=3.78e-4, FDR=0.019), cell-cell junction, synaptic membrane, and substance-related disorder gene sets. Separately, rare-variant burden flagged 55 genome-wide-significant genes (cell-cell junction, ECM, dynein/microtubule). **Did not cross-reference migraine GWAS or CGRP genes.**
**Novelty implication**: Most important paper for our framing — first systematic POTS genomics, and the estrogen-response enrichment is a direct bridge to estrogen-CGRP biology (see Krause 2021).

### `mueller_2022_potsmigraine`
Mueller & Robinson-Papp 2022, *Headache* 62(7):792–800. doi:10.1111/head.14365.
**Key claim**: Migraine is the most common comorbidity in POTS; proposed mechanistic overlaps are sympathetic dysregulation, hemodynamic alterations, and central sensitization. **Does not center CGRP as the bridging mechanism.**
**Novelty implication**: Strongest existing POTS-migraine review; CGRP omission is exactly our gap.

### `ray_2022_metaanalysis`
Ray et al. 2022, *Cephalalgia* 42(11-12):1274–1287. doi:10.1177/03331024221095153.
**Key claim**: Pooled meta-analytic prevalence of migraine in POTS = 36.8% (95% CI 2.9–70.7%) across 23 studies; high between-study heterogeneity.
**Novelty implication**: Best quantitative anchor for comorbidity claim; wide CI itself motivates better mechanistic understanding.

### `zhang_2021_hrv`
Zhang et al. 2021, *Front Neurol* 12:647092. doi:10.3389/fneur.2021.647092.
**Key claim**: Migraineurs show reduced HRV (lower SDNN, triangular index); parasympathetic dysfunction predominates ictally, sympathetic abnormalities also present.
**Novelty implication**: Autonomic dysfunction is a measurable feature of migraine itself, independent of comorbid POTS.

### `bendik_2011_hypermobility`
Bendik et al. 2011, *Cephalalgia* 31(5):603–613. doi:10.1177/0333102410392606.
**Key claim**: Migraine prevalence 75% in JHS/hEDS women vs 43% controls; earlier onset, more frequent.
**Novelty implication**: Quantitative anchor for hEDS-migraine link.

### `halverson_2023_heds`
Halverson et al. 2023, *Genet Med Open* 1(1):100812. doi:10.1016/j.gimo.2023.100812.
**Key claim**: hEDS co-diagnosis prevalence: anxiety 75%, depression 68%, **migraine 67%, POTS 61%**, IBS 57%; mean 10.45 co-diagnoses/individual.
**Novelty implication**: Largest dataset documenting the hEDS-POTS-migraine-MCAS phenotypic cluster.

### `krause_2021_estrogen` ⭐
Krause et al. 2021, *Nat Rev Neurol* 17(10):621–633. doi:10.1038/s41582-021-00544-2.
**Key claim**: Estrogen withdrawal disinhibits trigeminovascular CGRP expression and release. **Direct molecular bridge from sex-hormone biology to CGRP-driven migraine.**
**Novelty implication**: Half of the triangulation — connects female sex bias (shared by POTS and migraine) to CGRP. Combined with Qu 2025's "early estrogen response" enrichment in POTS rare variants, the estrogen→CGRP axis becomes the natural mechanistic bridge our paper articulates.

### Bucket C synthesis
POTS-migraine comorbidity is robust but quantitatively unstable: meta-analytic pooled migraine prevalence in POTS is 36.8% (CI 2.9–70.7%; Ray 2022); large clinic cohorts report 65–67% (Halverson 2023, Mueller 2022). Both share female predominance and a connective-tissue cluster (hEDS-POTS-migraine-MCAS; Halverson 2023, Bendik 2011), and migraineurs themselves have measurable autonomic dysregulation (Zhang 2021). POTS genetics was nearly absent until 2025: prior to Qu 2025, only one validated monogenic case (*SLC6A2*; Shannon 2000) and one candidate-gene signal (*GNB3*; Nakao 2012). **Qu 2025 is the first systematic POTS GWAS+WES — no common variants reached GWS, but over-representation analysis on nominally-significant common variants recovered HALLMARK_ESTROGEN_RESPONSE_EARLY (P=3.78e-4, FDR=0.019); rare-variant burden separately identified 55 genome-wide-significant genes.** This connects directly to Krause 2021's estrogen→CGRP mechanism. **No published autonomic-dysfunction-in-migraine review centers CGRP as the bridging mechanism** (Mueller 2022 does not, Krause 2021 does not address autonomic dysfunction).

---

## D. Computational pathway and network methods

### `deleeuw_2015_magma`
de Leeuw et al. 2015, *PLoS Comp Biol* 11(4):e1004219.
**Key claim**: MAGMA — gene-level + gene-set analysis from GWAS summary stats with multiple-regression aggregation accounting for LD; competitive gene-set tests with covariate adjustment.
**Novelty implication**: Canonical method for converting migraine GWAS summary stats into gene-level scores.

### `watanabe_2017_fuma`
Watanabe et al. 2017, *Nat Commun* 8:1826.
**Key claim**: FUMA — wraps MAGMA + positional/eQTL/chromatin SNP-to-gene mapping + GTEx tissue enrichment + MsigDB pathway tests in a reproducible pipeline.
**Novelty implication**: Turnkey pipeline; report version + parameters.

### `raudvere_2019_gprofiler`
Raudvere et al. 2019, *Nucleic Acids Res* 47(W1):W191–W198.
**Key claim**: g:Profiler suite for over-representation against GO/KEGG/Reactome/WikiPathways/TF/miRNA/HPA/HPO with hierarchical g:SCS correction tailored to nested ontologies.
**Novelty implication**: Preferred enrichment engine — g:SCS more honest than naive BH on redundant GO branches.

### `kuleshov_2016_enrichr`
Kuleshov et al. 2016, *Nucleic Acids Res* 44(W1):W90–W97.
**Key claim**: Enrichr — Fisher-style enrichment with combined-score ranking against >180 libraries (ontologies, pathways, drugs, TFs, cell types, disease signatures, ARCHS4, DSigDB, GTEx).
**Novelty implication**: Secondary check for drug-perturbation/tissue libraries g:Profiler does not carry.

### `guney_2016_proximity` ⭐
Guney et al. 2016, *Nat Commun* 7:10331.
**Key claim**: Defines and benchmarks network-proximity measures (closest, shortest, kernel, centre, separation) between drug-target sets and disease-gene sets in the human interactome with degree-preserving randomization null. "Closest" distance with degree-binned reference sets best discriminates effective from palliative drugs.
**Novelty implication**: Directly defines our migraine-vs-CGRP-pathway proximity test; published code (`toolbox`/`proximity`) is the de-facto standard.

### `menche_2015_diseaseome`
Menche et al. 2015, *Science* 347(6224):1257601.
**Key claim**: Mathematical conditions for disease-module identifiability in incomplete PPI networks; introduces network separation S_AB. Overlapping modules (S_AB < 0) predict comorbidity, co-expression, and symptom similarity.
**Novelty implication**: Report S_AB alongside proximity. Identifiability caveat — needs ~25–50% of true module genes — goes in Limitations.

### `maslov_2002_rewiring`
Maslov & Sneppen 2002, *Science* 296(5569):910–913.
**Key claim**: Degree-preserving edge-swap rewiring as the canonical null for PPI analyses.
**Novelty implication**: Gold-standard null; Guney 2016 degree-binned random gene sets are the computationally cheaper equivalent.

### `szklarczyk_2023_string`
Szklarczyk et al. 2023, *Nucleic Acids Res* 51(D1):D638–D646.
**Key claim**: STRING v12 — physical + functional associations across ~12K organisms; per-edge confidence; integrates text-mining, experiments, databases, co-expression, genomic context.
**Novelty implication**: Includes text-mining edges → **literature bias**. Use physical-only or high-confidence subnetwork; report as sensitivity.

### `oughtred_2021_biogrid`
Oughtred et al. 2021, *Protein Sci* 30(1):187–200.
**Key claim**: BioGRID — ~1.93M curated physical/genetic/chemical interactions; PubMed-anchored evidence codes; low- vs high-throughput stratification.
**Novelty implication**: Strong candidate for **physical** PPI backbone (no text-mining inflation).

### `deltoro_2022_intact`
Del Toro et al. 2022, *Nucleic Acids Res* 50(D1):D648–D653.
**Key claim**: IntAct (IMEx consortium) — >1M binary interactions with PSI-MI-curated detail (method, role, kinetics).
**Novelty implication**: Second curated physical-only source; union with BioGRID for sensitivity analysis.

### `gtex_2020_atlas`
GTEx Consortium 2020, *Science* 369(6509):1318–1330.
**Key claim**: GTEx v8 — cis/trans-eQTL + sQTL maps from 15,201 RNA-seq samples / 49 tissues / 838 donors.
**Novelty implication**: eQTL substrate inside FUMA; tissue-enrichment source. Lacks trigeminal ganglion / autonomic ganglion (proxy with brain stem, tibial nerve).

### `goeman_2007_genesets`
Goeman & Bühlmann 2007, *Bioinformatics* 23(8):980–987.
**Key claim**: Distinguishes "competitive" vs "self-contained" null hypotheses; gene-sampling competitive tests assume gene-gene independence (almost always violated).
**Novelty implication**: When over-representing CGRP-pathway among migraine GWAS hits, prefer MAGMA's regression-based competitive test (LD-adjusted) over naive hypergeometric.

### Bucket D synthesis
**Recommended stack**: (1) MAGMA via FUMA on published migraine GWAS summary stats, with GTEx v8 eQTL SNP-to-gene mapping and tissue stratification → versioned reproducible migraine gene list. (2) g:Profiler primary (g:SCS); Enrichr secondary for drug/tissue libraries; flag Goeman & Bühlmann caveat — MAGMA gene-set test is statistically preferred primary. (3) Guney "closest" proximity d_c between migraine gene set and CGRP-pathway gene set, with Menche separation S_AB reported alongside.

**PPI database choice**: Default to **BioGRID ∪ IntAct (physical only)**; STRING (high-confidence physical, text-mining disabled, score ≥700) as sensitivity. STRING text-mining channel artificially shortens distance between any two well-studied gene sets — and CGRP genes are heavily literature-studied — so it would bias toward false-positive proximity.

**Null model**: Degree-preserving randomization — 1000 size-matched random gene sets where each random gene is drawn from same degree bin as its real counterpart (Guney 2016 protocol; cheaper equivalent of Maslov–Sneppen edge rewiring). Pure label permutation is **not** degree-preserving and inflates significance for high-degree disease genes. Also include matched-control diseases (psoriasis, schizophrenia) to demonstrate specificity.

**Limitations to flag**: interactome incompleteness (Menche identifiability threshold); ascertainment bias (well-studied genes accumulate edges); competitive-test independence violation (Goeman 2007); GTEx tissue gaps (no trigeminal/autonomic ganglion samples); reproducibility (record FUMA/MAGMA/g:Profiler/Enrichr/STRING/BioGRID/IntAct versions + access dates + random seeds).

---

## E. Existing triangulation work — gap analysis

### `siewert_2020_crosstrait`
Siewert et al. 2020, *Int J Epidemiol* 49(3):1022–1031. doi:10.1093/ije/dyaa050.
**Key claim**: LDSC genetic correlations between migraine and 47 UK Biobank phenotypes + external GWAS. Found rg=0.10 with diastolic BP (P=5.4e-5), rg=0.063 with systolic BP, plus heart disease, lipids, platelet/WBC counts. **Did NOT include HRV, orthostatic intolerance, syncope, or POTS-related phenotypes.** No PPI proximity, no CGRP-pathway gene-set test.
**Overlap with our plan**: PARTIAL — establishes migraine-vascular genetic correlation by LDSC, but stops at BP and never reaches autonomic/HRV/POTS or CGRP PPI proximity.

### `guo_2020_bp_migraine`
Guo et al. 2020, *Nat Commun* 11:3368. doi:10.1038/s41467-020-17002-0.
**Key claim**: Cross-phenotype meta-analysis of migraine and BP GWAS; 14 shared loci. **ADRA2B (alpha-2B adrenergic receptor) is the only adrenergic-receptor gene called out.** No formal autonomic gene-set enrichment, no CGRP PPI proximity, no POTS framing.
**Overlap with our plan**: ADJACENT — surfaces a single adrenergic-receptor shared locus en route to BP cross-phenotype, no permutation test against an adrenergic gene set.

### `zhang_2024_druggable`
Zhang et al. 2024, *J Headache Pain* 25:100. doi:10.1186/s10194-024-01805-3.
**Key claim**: Druggable-genome MR + colocalization on migraine GWAS; 21 druggable migraine-associated genes. PPI used as downstream visualization. CALCB/CALCRL/RAMP1/RAMP3 referenced only as comparator known drug targets, not as a tested complex. **No autonomic gene-set test, no PPI proximity to CGRP-complex anchor, no POTS framing.**
**Overlap with our plan**: ADJACENT — uses GWAS+PPI machinery for target discovery, not gap-test design.

### `lin_2024_syncope_migraine`
Lin et al. 2024, *J Clin Neurol* 20(6):599–609. doi:10.3988/jcn.2024.0156.
**Key claim**: GWAS of syncope within migraine subtypes; tested only three pre-selected variants (TRPM8, GFRA1, LRP1), all non-significant.
**Overlap with our plan**: TANGENTIAL — touches migraine + autonomic-symptom genetics but only via three candidate SNPs.

### `rappoport_2024_sns_migraine` ⭐ (closest conceptual cousin)
Rappoport 2024, arXiv:2408.06780.
**Key claim**: Verbal/mechanistic theory linking SNS partial saturation, baroreceptor hyperexcitability, and downstream CGRP release. Discusses TRPM8, PRDM16, ATP1A2 candidate genes from prior literature. **Performs no GWAS overlap, no gene-set enrichment, no PPI proximity analysis.**
**Overlap with our plan**: ADJACENT (conceptual) — proposes the same autonomic-CGRP linkage we test, but does not perform the computational test. Cite as motivating prior art.

### `liu_2025_precision_target`
Liu X et al. 2025, *Molecules* 30(19):3921. doi:10.3390/molecules30193921.
**Key claim**: SMR/HEIDI + colocalization + PPI + enrichment pipeline on migraine GWAS; prioritized 8 genes; computationally proposed 41 repurposable drugs. PPI used to score interaction with known CGRP-/5-HT-pathway drug targets, but not as a network-proximity hypothesis test. No autonomic gene-set test; no POTS framing.
**Overlap with our plan**: ADJACENT — uses PPI in a related way (interaction with CGRP-pathway drug targets) but the construct is target prioritization, not a permutation-controlled set-overlap or proximity-vs-null test.

### Negative search evidence

The following targeted searches returned **no paper**:
- migraine GWAS × adrenergic/sympathetic/autonomic gene set with permutation null
- migraine GWAS hits closer to CALCA/CALCB/CALCRL/RAMP1 in PPI than chance
- migraine × POTS LDSC genetic correlation (Qu 2025 has not been cross-trait-analyzed with migraine)
- migraine + heart-rate-variability LDSC
- trigeminovascular network proximity to autonomic gene set

---

## Gap-analysis verdict

**PARTIALLY NOVEL — leaning strongly novel.**

Each leg of the proposed triangle exists in the literature individually, but no single computational object combines them, and the specific tests we plan have not been published:

1. **Migraine GWAS + generic pathway enrichment** is fully covered (Hautakangas 2022 Supp Tables 15–16 via MAGMA/DEPICT; Choquet 2021 via FUMA/DEPICT). **Neither tested an autonomic / adrenergic-receptor / sympathetic-nervous-system gene set as a focused hypothesis with a permutation null distinct from generic MSigDB enrichment.**

2. **Migraine + cardiovascular genetic correlation via LDSC** is covered (Siewert 2020 for BP/heart-disease; Guo 2020 for BP cross-phenotype with ADRA2B locus surfacing). **No published LDSC genetic correlation between migraine and an autonomic / HRV / orthostatic / POTS phenotype exists** — Siewert's 47-trait panel stops at BP, and the only existing POTS GWAS (Qu 2025) is too underpowered for genome-wide significance and has not been cross-trait-analyzed with migraine.

3. **Migraine drug-target / PPI / network-pharmacology** is partly covered (Zhang 2024, Liu X 2025) but every published instance is shaped as **target prioritization**, not as a hypothesis test of "are migraine GWAS hits closer to the CGRP receptor complex (CALCA/CALCB/CALCRL/RAMP1) than degree-matched random gene sets?" with a Guney-style permutation null.

4. **POTS-comorbidity-motivated framing**: nothing exists. Rappoport 2024 (arXiv) proposes the SNS–CGRP linkage verbally but performs no computation; no paper crosses migraine GWAS, an autonomic gene set, and CGRP PPI proximity in a single object.

5. **Bonus convergence**: Qu 2025 POTS GWAS over-representation analysis (on nominally-significant common variants) recovered HALLMARK_ESTROGEN_RESPONSE_EARLY (FDR=0.019), which (via Krause 2021) is mechanistically continuous with estrogen-modulated CGRP biology. This convergence has not been articulated in any single paper.

### What our paper would add
- A focused permutation-null test of migraine GWAS-prioritized genes against curated autonomic / adrenergic gene sets (rather than scanning thousands of MSigDB sets).
- A network-proximity test (Guney 2016 statistic + degree-preserving null) of migraine GWAS gene set vs the CGRP receptor complex (CALCRL, RAMP1, RAMP2, RAMP3, CALCA, CALCB) as the disease-anchor module.
- POTS-comorbidity framing as the motivating clinical observation, with the estrogen→CGRP axis (Qu 2025 + Krause 2021) as the natural mechanistic bridge.

### Recommendation
**Proceed with the project.** Cite Hautakangas 2022, Choquet 2021, Siewert 2020, Guo 2020, Zhang 2024, and Rappoport 2024 in the introduction to make the gap explicit. Note specifically that (a) MAGMA/DEPICT enrichments did not interrogate autonomic gene sets, (b) no LDSC has been computed for migraine vs HRV/orthostatic/POTS phenotypes, and (c) prior PPI/network-pharmacology work on migraine has been target-prioritization, not CGRP-anchored proximity hypothesis testing.
