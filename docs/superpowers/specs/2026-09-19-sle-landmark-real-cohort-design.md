# SLE Landmark 真实队列比较设计

## 目标

在固定 STRING 高置信网络、Reactome、RWR、SVD 和 SLE 差异分析结果的前提下，只改变 Landmark 选择方法，比较 `degree_high`、`weighted_betweenness_high` 和 `weighted_kmedoids_high` 在三个独立 SLE 队列中的下游复现性。

## 文件结构

实现保持为用户要求的三个 R 脚本：

1. `scripts/real_cohort/01_run_sle_landmark_benchmark.R`：唯一运行入口。集中声明参数和路径，加载已有 Landmark，构建三套 fit，依次运行 Weighted-GSEA 与 Concise-ORA，并调用已有跨队列比较函数。
2. `scripts/real_cohort/02_sle_landmark_helpers.R`：只补现有代码缺少的胶水函数，包括用固定 Landmark 构建 `gsemb_embedding`、整理 ranked DEG、规范化富集结果、组织结果映射和汇总两条路线。
3. `scripts/real_cohort/03_sle_landmark_beginner_walkthrough.R`：可阅读的中文教学脚本，按运行顺序解释输入、参数、中间对象、输出和结论边界；默认不启动昂贵计算。

函数测试放在 `tests/testthat/test-sle-landmark-real-cohort.R`。

## 复用边界

- 直接复用包内 `gsemb_compute_node_landmark_features()`、`gsemb_fit_gene_embedding()`、`gsemb_fit_set_gaussians_from_members()`，不重写 RWR、SVD 或高斯拟合。
- 直接复用 `gsemb_weighted_gsea()` 与 `gsemb_concise_gene_sets()`。
- 直接复用 `step3/step3_modular/01_ora_from_deg.R` 和 `02_compare_ora.R` 中的 ORA、Jaccard、overlap coefficient 与 Spearman 实现。
- 直接读取 worktree 已有 degree、betweenness、K-Medoids Landmark RDS；不重新运行耗时选点。
- 不新增依赖，不修改 `gsemb_fit()` 公共 API。

## 两条分析路线

### 路线 A：Weighted-GSEA（主分析）

每个 Landmark 方法生成一套 fit。每个 SLE 队列以命名 `log2FoldChange` 向量作为 ranked statistic，调用现有 `gsemb_weighted_gsea()`。结果统一为含 `ID`、`pvalue`、`p.adjust` 的表，再调用现有跨队列比较函数计算完整排名 Spearman 和 Top-10/20/50 Jaccard。

### 路线 B：Concise-ORA（敏感性分析）

每套 fit 使用相同的原始 Reactome 成员，通过固定的 Concise 参数生成方法特异的 Concise 通路。随后对三个队列使用相同 DEG 阈值运行现有 ORA，再调用同一套跨队列比较函数。该路线同时包含“Landmark 改变嵌入”和“嵌入改变 Concise 成员”两段效应，因此只作为敏感性分析，不取代主分析。

## 固定比较条件

- PPI：STRING v12，`combined_score > 700`。
- Landmark 数：`k = 128`。
- 方法：`degree_high`、`weighted_betweenness_high`、`weighted_kmedoids_high`。
- 队列：`GSE50772`、`GSE61635`、`GSE72509`。
- Reactome、SVD 维数、RWR 参数、随机种子、富集参数在三种方法间完全相同。
- 主终点：三队列两两完整通路排名 Spearman 的方法内均值。
- 辅助终点：Top-10/20/50 Jaccard 与 overlap coefficient。
- 显著通路数量只作描述，不作为选择赢家的主标准。

## 输出

输出根目录按时间戳创建，包含：

- `fits/`：三套 fit RDS；
- `weighted_gsea/`：9 张队列×方法结果表与跨队列汇总；
- `concise_ora/`：三套 Concise 集合、9 张 ORA 表与跨队列汇总；
- `run_meta.txt`：代码路径、输入路径、参数和随机种子；
- `method_summary.csv`：按路线和方法汇总的复现指标。

昂贵结果若已存在则默认复用，显式 `overwrite <- TRUE` 才覆盖。

## 结论边界

本实验可以支持“某种 Landmark 表示在当前三个 SLE 队列中产生更高的下游复现性”。它不能证明 Landmark 节点是 SLE 金标准，也不能把三个队列上的结果外推到所有疾病。Concise-ORA 的差异不能单独归因于 Landmark，因为通路成员也随 fit 改变。
