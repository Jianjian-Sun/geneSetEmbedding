# Landmark 真实队列基准：SLE 与 BRCA 结果整理

- **日期**：2026-09-20  
- **工作树**：`geneSetEmbedding/.worktrees/weighted-degree-landmarks`  
- **状态**：本轮已重跑；下文数字均来自当次产物表，非申请书预期  

| 疾病 | 输出目录 | 入口脚本 |
|------|----------|----------|
| SLE | `test_results/real_cohort/sle_landmark_benchmark_v1/` | `scripts/real_cohort/01_run_sle_landmark_benchmark.R` |
| BRCA | `test_results/real_cohort/brca_landmark_benchmark_v1/` | `scripts/real_cohort/04_run_brca_landmark_benchmark.R` |

总表：

- SLE：`sle_landmark_benchmark_v1/method_summary.csv`
- BRCA：`brca_landmark_benchmark_v1/method_summary.csv`

---

## 1. 实验在比什么

在**固定** STRING 高置信网络（`combined_score > 700`）、Reactome、RWR、SVD、富集参数的前提下，只换 Landmark 选法，比较三种代表：

| 方法名 | 含义 |
|--------|------|
| `weighted_degree_high` | 高置信图上按加权度（strength）Top-128 |
| `weighted_betweenness_high` | 加权介数 Top-128（读既有 RDS） |
| `weighted_kmedoids_high` | 加权 K-Medoids 128 中心（读既有 RDS） |

共同设定（两病一致）：`k=128`，`embedding_dim=64`，`alpha=0.5`，`seed=20260919`，Weighted-GSEA `nperm=1000`。

Landmark / 基因嵌入与疾病无关。BRCA 跑次复用了 SLE 目录下已缓存的 fit 与 Concise 基因集，只更换各队列 DEG / logFC 做下游富集。

### 两条路线

1. **主分析 Weighted-GSEA**  
   输入为全基因组连续 `log2FoldChange`（**不**先按 padj 筛 DEG）。Landmark → 嵌入 → 通路软权重 → 加权富集。  
   **主终点**：三队列两两通路完整排名的 Spearman，再对方法取均值（`mean_rank_spearman`）。

2. **敏感性 Concise-ORA**  
   嵌入后对 Reactome 做 Concise，再用硬阈值显著 DEG 做普通 ORA。  
   同时含“Landmark 改嵌入”和“嵌入改 Concise 成员”，**不能单独归因于 Landmark**。

### 指标含义

| 指标 | 含义 | 注意 |
|------|------|------|
| `mean_rank_spearman` | 跨队列通路**完整排名**一致性 | 主终点；不依赖 FDR 截断 |
| `mean_jaccard` | 跨队列 **FDR&lt;0.05 显著通路集合**重叠 | 一侧显著数≈0 时会塌成 0 / NA |
| `mean_overlap_coefficient` | 交集 / min(\|A\|,\|B\|) | 比 Jaccard 更宽松 |

这些指标量的是**跨队列复现性**，不是“生物学金标准正确性”，也不能证明 Landmark 节点本身是疾病调控基因。

---

## 2. 队列与输入同向性（背景）

### SLE

- 队列：`GSE50772`，`GSE61635`，`GSE72509`  
- DEG：`E:/genesetembedding/step3/SLE/SLE_results/*_deg_all.csv`  
- 平台差异大：两芯片（GPL570）+ 一 RNA-seq（GPL16791）；PBMC / blood / whole blood  

全基因组 logFC Spearman（不筛 padj）：

| 队列对 | Spearman |
|--------|----------|
| GSE50772 vs GSE61635 | 0.177 |
| GSE50772 vs GSE72509 | 0.241 |
| GSE61635 vs GSE72509 | 0.437 |
| **均值** | **≈ 0.29** |

硬阈值显著 DEG Jaccard 均值约 **0.07**。输入本身跨队列较散，下游复现天花板偏低。

### BRCA

- 队列：`GSE109169`，`GSE29044`，`TCGA_BRCA`  
- DEG：`E:/genesetembedding/step3/BRCA/BRCA_results/*_deg_all.csv`  

全基因组 logFC Spearman：

| 队列对 | Spearman |
|--------|----------|
| GSE109169 vs GSE29044 | 0.751 |
| GSE109169 vs TCGA_BRCA | 0.688 |
| GSE29044 vs TCGA_BRCA | 0.661 |
| **均值** | **≈ 0.70** |

输入明显比 SLE 同向；更适合作为“跨队列复现”考场。

---

## 3. SLE 结果

来源：`sle_landmark_benchmark_v1/method_summary.csv`（跑次时间见该目录 `run_meta.txt`：2026-09-20 16:22:18）。

### 3.1 主分析 Weighted-GSEA

| 方法 | mean_rank_spearman | mean_jaccard (all_sig) |
|------|--------------------|-------------------------|
| **weighted_kmedoids_high** | **0.332** | 0.034 |
| weighted_degree_high | 0.136 | 0 |
| weighted_betweenness_high | 0.110 | 0 |

**说明：**

- 按预设主终点，SLE 上 **kmedoids 相对最好**，但绝对值不高（约 0.33）。  
- degree / betweenness 的 Jaccard≈0：多数队列对上 FDR 显著通路数为 0（例如 degree 在 GSE50772、GSE72509 均为 0 显著，仅 GSE61635 大量显著）。显著集合重叠指标在本轮**几乎失效**。  
- 完整排名 Spearman 仍可计算；kmedoids 三对中 GSE50772–GSE72509 较高（约 0.52），拉高了均值。

### 3.2 敏感性 Concise-ORA

| 方法 | mean_rank_spearman | mean_jaccard (all_sig) |
|------|--------------------|-------------------------|
| **weighted_kmedoids_high** | **0.401** | **0.467** |
| weighted_betweenness_high | 0.315 | 0.316 |
| weighted_degree_high | 0.308 | 0.213 |

方向与主分析一致（kmedoids 最好）。Jaccard 0.47 表示三对显著通路集合平均约一半重叠（例如 8∩12→7 等），是集合重叠意义上的相对优势，**不是**“已证明稳定”，且含 Concise 混淆。

### 3.3 SLE 可写 / 不可写

**可写：** 在本设定下，`weighted_kmedoids_high` 使三个 SLE 队列的下游通路排名（及 Concise-ORA 显著集合）**相对**另两种更可复现。  

**不可写：** Landmark 是 SLE 金标准；绝对值已“够稳”；已证明更好保留生物学通路；GSEA Jaccard 失败等于方法完全无效。

---

## 4. BRCA 结果

来源：`brca_landmark_benchmark_v1/method_summary.csv`（`run_meta.txt`：2026-09-20 18:05:07）。

### 4.1 主分析 Weighted-GSEA

| 方法 | mean_rank_spearman | mean_jaccard (all_sig) |
|------|--------------------|-------------------------|
| **weighted_betweenness_high** | **0.736** | **0.343** |
| weighted_kmedoids_high | 0.710 | 0.161 |
| weighted_degree_high | 0.432 | 0 |

**说明：**

- BRCA 上 Spearman 整体远高于 SLE（约 0.43–0.74），与输入 logFC 更同向一致。  
- 主终点赢家变为 **betweenness**（略高于 kmedoids）。  
- degree 的显著通路 Jaccard 仍为 0（多队列显著数失衡）；但其 Spearman（0.43）仍明显低于另两种。  
- betweenness 三对 Spearman 约 0.68–0.78，且三队列均有一定数量显著通路，Jaccard 也可读。

### 4.2 敏感性 Concise-ORA

| 方法 | mean_rank_spearman | mean_jaccard (all_sig) |
|------|--------------------|-------------------------|
| weighted_betweenness_high | 0.558 | 0.560 |
| weighted_degree_high | 0.554 | 0.566 |
| weighted_kmedoids_high | 0.542 | 0.565 |

三种方法非常接近；显著通路数量也远多于 SLE（每队列约 110–150 条）。**不能**用 Concise-ORA 在 BRCA 上强行区分 Landmark 优劣。

### 4.3 BRCA 可写 / 不可写

**可写：** 在输入更同向的 BRCA 三队列上，Weighted-GSEA 跨队列排名一致性整体提高；`weighted_betweenness_high` 与 `weighted_kmedoids_high` 明显优于 `weighted_degree_high`，二者接近、betweenness 略高。  

**不可写：** 已选出全局最优 Landmark；BRCA 结论自动外推到 SLE 或所有疾病。

---

## 5. 两病对照与综合解读

| 项目 | SLE | BRCA |
|------|-----|------|
| 输入 logFC Spearman 均值 | ≈ 0.29 | ≈ 0.70 |
| GSEA 主终点最优 | kmedoids (0.33) | betweenness (0.74) |
| GSEA 第二 | degree (0.14) | kmedoids (0.71) |
| GSEA 最弱 | betweenness (0.11) | degree (0.43) |
| Concise-ORA | kmedoids 明显更好 | 三种方法几乎持平 |
| GSEA 显著 Jaccard | 基本失效 | betweenness 可读；degree 仍塌 |

### 综合判断（证据等级：单次完整跑通，E2–E3）

1. **跨队列复现高度依赖疾病队列同质性。** SLE 输入散 → 所有方法 Spearman 偏低；BRCA 输入齐 → 天花板升高。  
2. **没有单一 Landmark 在两病主终点上同时“完胜”。** SLE 偏 kmedoids；BRCA 偏 betweenness（与 kmedoids 接近）。  
3. **`weighted_degree_high` 在两病 Weighted-GSEA 主终点上均最弱或明显偏弱。**  
4. **显著通路 Jaccard 不宜作为本轮主判据**（尤其 degree；SLE 几乎全军覆没）。应以完整排名 Spearman 为主，Jaccard 仅作描述。  
5. **当前指标只能做相对复现比较**，不能证明“保留真实生物学通路”。若要正确性，需 curated 通路/基因回收或随机 Landmark 基线等另设标准。

---

## 6. 产物结构（两目录相同）

```text
*_landmark_benchmark_v1/
  method_summary.csv          # 最先看的总表
  run_meta.txt                # 参数与 provenance
  fits/                       # 三套 gsemb_embedding RDS
  weighted_gsea/              # 9 张队列×方法结果 + cross_cohort/
  concise_ora/                # Concise RDS、ORA 结果与图 + cross_cohort/
```

跨队列明细 CSV 前缀：

- SLE：`SLE_landmark_weighted_gsea_*`，`SLE_landmark_concise_ora_*`
- BRCA：`BRCA_landmark_weighted_gsea_*`，`BRCA_landmark_concise_ora_*`

---

## 7. 一句话结论

- **SLE**：考场难；相对赢家是 **kmedoids**，但 Spearman 绝对值有限。  
- **BRCA**：考场更齐；**betweenness ≈ kmedoids ≫ degree**（主看 Weighted-GSEA Spearman）。  
- **共同**：degree 作为默认锚点在真实队列复现上不占优；最终默认选哪套，还需预注册规则 +（建议）随机基线 / 第二疾病复核，不宜单凭一轮相对排序定案。
