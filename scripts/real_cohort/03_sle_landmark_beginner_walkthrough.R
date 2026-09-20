# =============================================================================
# 给初学者看的说明：三种 Landmark 如何在真实 SLE 队列中比较
# =============================================================================
#
# 这个文件的主要用途是“边读边理解”，不是重新实现分析算法。
# 真正运行分析的入口是：
#   scripts/real_cohort/01_run_sle_landmark_benchmark.R
#
# 为防止你只是阅读本文件时误启动耗时计算，默认设为 FALSE。
RUN_ANALYSIS <- FALSE

# =============================================================================
# 一、这次比较到底在比较什么？
# =============================================================================
#
# 先在每个方法家族内部选一个代表：
#
#   Weighted degree -> weighted_degree_high
#   Betweenness     -> weighted_betweenness_high
#   K-Medoids       -> weighted_kmedoids_high
#
# 三者都使用同一张 STRING 高置信网络（combined_score > 700），Landmark
# 数量也都固定为 128。这样最终差异才可以主要归因于 Landmark 方法。
#
# Landmark 不是“真正的 SLE 致病基因名单”，也不是要直接做 ORA 的基因集。
# 它们是网络坐标锚点。每个普通基因都会得到“它相对于 128 个 Landmark
# 的扩散位置”，随后这些位置再经过 SVD 变成基因嵌入。
#
# 简化的数据流：
#
#   STRING PPI
#      -> 选 128 个 Landmark
#      -> 从每个 Landmark 做 RWR
#      -> 每个基因得到 128 维扩散特征
#      -> SVD 得到基因嵌入
#      -> Reactome 通路得到高斯表示
#      -> 用 SLE 的差异表达结果做富集
#

# =============================================================================
# 二、为什么有两条路线？
# =============================================================================
#
# 路线 A：Weighted-GSEA（主分析）
# --------------------------------
# 它把每个 SLE 队列中所有可用基因的 log2FoldChange 作为连续排序统计量，
# 再利用每套 Landmark fit 产生的“基因到通路软权重”计算富集。
#
# Landmark 变化会改变 RWR 特征和基因嵌入，进而改变软权重和最终富集。
# 所以这条路线最直接地回答：
#
#   只改变 Landmark，哪个方法产生的 SLE 通路结果更容易跨队列复现？
#
# 路线 B：Concise-ORA（敏感性分析）
# ----------------------------------
# 每套 Landmark fit 先从每条原始 Reactome 通路的成员中选一个更精简的
# Concise 子集，再用同一批 SLE DEG 做普通 ORA。
#
# 注意：这里不是“把 Landmark 再做 Concise”。Concise 筛的是 Reactome
# 通路成员，Landmark 只是间接改变了基因嵌入和筛选次序。
#
# 这条路线同时包含两段影响：
#
#   Landmark 改变嵌入 + 嵌入改变 Concise 通路成员
#
# 因此它只能作为敏感性分析，不能替代主分析，也不能把所有差异都说成
# Landmark 单独造成的。
#

# =============================================================================
# 三、三个脚本分别做什么？
# =============================================================================
#
# 01_run_sle_landmark_benchmark.R
#   - 你真正需要运行的主脚本。
#   - 顶部集中列出所有参数和路径。
#   - 加载已有 Landmark、构建或复用 fit、运行两条路线、写出结果。
#
# 02_sle_landmark_helpers.R
#   - 主脚本自动调用，不需要你单独运行。
#   - 只补包里原来缺少的连接步骤，例如“使用指定 Landmark 构造 fit”。
#   - RWR、SVD、高斯拟合、GSEA、ORA、Jaccard 都继续调用已有实现。
#
# 03_sle_landmark_beginner_walkthrough.R
#   - 就是你正在读的文件。
#   - 解释为什么这么做、每个结果怎么看、什么结论不能扩大。
#

# =============================================================================
# 四、主脚本使用了哪些现有数据？
# =============================================================================
#
# 三个队列：
#   E:/genesetembedding/step3/SLE/SLE_results/GSE50772_deg_all.csv
#   E:/genesetembedding/step3/SLE/SLE_results/GSE61635_deg_all.csv
#   E:/genesetembedding/step3/SLE/SLE_results/GSE72509_deg_all.csv
#
# 每张表至少使用：
#   gene_symbol     基因名称
#   log2FoldChange  病例相对对照的连续效应值，供 Weighted-GSEA 排序
#   padj            ORA 路线筛 DEG 时使用的多重校正 P 值
#
# Reactome：
#   E:/genesetembedding/step3/reactome_gene_sets.RData
#
# 已有 Landmark：
#   Weighted degree 因旧脚本没保存专用 RDS，直接在已缓存的 high 图上快速重算。
#   Betweenness 直接读取 betweenness_landmarks.rds。
#   K-Medoids 直接读取 kmedoids_landmarks.rds。
#
# 几何缓存图与 Landmark 使用 STRING protein ID（9606.ENSP...）。
# 主脚本会用 step1/STRING/9606.protein.info.v12.0.txt 映射到 gene symbol，
# 以便与 DEG / Reactome 对齐。这避免了重复运行耗时的 Betweenness 和 K-Medoids。
#

# =============================================================================
# 五、哪些参数可以改？
# =============================================================================
#
# 第一次正式比较建议保持主脚本默认值。尤其不要给三个方法设置不同参数。
#
# 可以统一修改：
#   gsea_nperm       置换次数。调试可先用 50；正式结果建议至少 1000。
#   top_k_values     比较 Top-10、20、50，可统一增减。
#   run_label        新实验使用新名称，避免与旧结果混在一起。
#   overwrite        FALSE 会复用结果；TRUE 会覆盖并重算。
#
# 暂时不要随意修改：
#   k                三种方法必须相同，默认 128。
#   embedding_dim    三种方法必须相同，默认 64。
#   alpha/tol/max_iter
#                    RWR 参数必须相同。
#   ppi_score_cutoff 三种方法必须使用同一 high 图。
#   concise 参数     Concise-ORA 中三种方法必须一致。
#

# =============================================================================
# 六、输出目录怎么看？
# =============================================================================
#
# 默认根目录：
#   test_results/real_cohort/sle_landmark_benchmark_v1/
#
# fits/
#   三套基因嵌入 fit。第二次运行会直接复用，省去 RWR/SVD 时间。
#
# weighted_gsea/
#   每个队列 × 每种方法一张结果表，共 9 张。
#   cross_cohort/ 中是三队列两两比较及汇总。
#
# concise_ora/
#   三套方法特异的 Concise 基因集，以及 9 张 ORA 结果表。
#   cross_cohort/ 中是三队列两两比较及汇总。
#
# method_summary.csv
#   最先看的总表。主要列：
#
#   route
#     weighted_gsea 是主分析；concise_ora 是敏感性分析。
#
#   method
#     三种 Landmark 代表方法。
#
#   pathway_mode
#     top10/top20/top50/all_sig。
#
#   mean_rank_spearman
#     三个队列对的完整通路排名相关均值。主分析优先看这个指标。
#
#   mean_jaccard
#     三个队列对的 Top-K 通路集合 Jaccard 均值。
#
#   mean_overlap_coefficient
#     交集除以两个集合中较小者的大小，对集合大小差异比 Jaccard 宽松。
#
# run_meta.txt
#   保存本次分析所有关键参数。汇报结果时应把它和结果表一起保留。
#

# =============================================================================
# 七、怎样判断哪种方法更好？
# =============================================================================
#
# 预先使用以下顺序，避免结果出来后挑对自己有利的指标：
#
# 1. 主看 weighted_gsea 路线的 mean_rank_spearman。
# 2. 再看 Top-20 和 Top-50 的 mean_jaccard 是否同方向支持。
# 3. 检查优势是否由三个队列对共同贡献，而不是只有一对特别高。
# 4. 看 concise_ora 是否提供方向一致的敏感性证据。
# 5. 显著通路数量只能描述，不能当作“越多越好”。
#
# 如果某方法只在 Concise-ORA 好、Weighted-GSEA 不好，应写成：
#   “优势依赖于 Concise 通路成员筛选，不能归因于 Landmark 本身。”
#
# 如果某方法在 Weighted-GSEA 的完整排名、Top-20 和 Top-50 都较高，才可写：
#   “该 Landmark 表示在当前三个 SLE 队列中保留了更可重复的疾病相关
#    通路信号。”
#

# =============================================================================
# 八、哪些结论不能说？
# =============================================================================
#
# 不能说：
#   - 这 128 个节点是 SLE 的真实 Landmark 金标准。
#   - Jaccard 高就证明所有基因都被完整保留。
#   - 在三个 SLE 队列最好，就在所有疾病都最好。
#   - 显著通路更多，所以方法更准确。
#
# 可以说：
#   - 在固定网络、参数和通路库条件下，某种 Landmark 方法使三个独立
#     SLE 队列的下游通路排序或 Top-K 通路集合具有更高复现性。
#

# =============================================================================
# 九、如何运行？
# =============================================================================
#
# 方式 1：打开主脚本，检查顶部路径和参数，然后在 R 中 source：
#
#   source("scripts/real_cohort/01_run_sle_landmark_benchmark.R")
#
# 方式 2：如果你确认要从这个教学文件启动，把下面的 FALSE 改成 TRUE，
# 然后 source 本文件。默认 FALSE 是为了避免误运行数小时任务。

if (isTRUE(RUN_ANALYSIS)) {
  source("scripts/real_cohort/01_run_sle_landmark_benchmark.R")
}
