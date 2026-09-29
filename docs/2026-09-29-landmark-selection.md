# Landmark 选择：三种加权算法与边表清洗

日期：2026-09-29。当前主分支接口以 `gsemb_select_landmarks()` / `gsemb_fit()` 为准。

Landmark 是扩散特征的参考节点。包从 PPI 邻接矩阵里选出 `k` 个节点，再对每个 landmark 做随机游走重启，得到基因嵌入。可选方法有四种：`degree`、`betweenness`、`kmedoids`、`random`。前三种都使用边权；`random` 是均匀随机对照，不看网络结构。

介数和 K-Medoids 的边代价定为 `1 / weight`。权重大的边更短。非零边权必须是有限正数；可以是 0–1、STRING 的 0–1000 分数或其他正数范围，无需强制归一化。

---

## 1. 三种算法

入口是 `gsemb_select_landmarks()`。`gsemb_fit()` 用同名参数 `landmark_method` 把选择交给这个函数，默认 `degree`。

### 1.1 `degree`：节点强度

对每个节点，把与它相连的边权加起来（邻接矩阵的行和），从大到小取前 `k` 个。

这是加权度，不是“有几条边”。一条权重 100 的边比十条权重 1 的边贡献更大。无向图在 `gsemb_build_graph()` 里已经做成对称矩阵，行和等于该节点的 strength。

并列时 `order()` 稳定排序，先出现在节点名顺序里的节点排在前面。

不使用 `seed`、`betweenness_cutoff`、`kmedoids_m`、`kmedoids_max_iter`。

### 1.2 `betweenness`：加权介数

一条边的长度是 `1 / weight`。节点的介数是：在所有节点对的最短路里，有多少条经过它。高介数节点更像网络上的桥。

实现调用 `igraph::betweenness()`，无向，边权为上述长度。算完后同样按分数从高到低取前 `k` 个。

在 STRING `combined_score > 700` 的人类图上，边长大约是 `1/700` 到 `1/999`。此前缓存图（重复方向的权重被相加）的所有有限最短路代价上界约 `0.015`，所以当时的 `cutoff=4` 没有截断任何路径。当前默认 `-1` 直接计算完整介数；新图的具体路径长度需要按该图复核。

### 1.3 `kmedoids`：加权 K-Medoids

目标是在图上放 `k` 个中心，使节点到最近中心的距离之和尽量小。距离同样是边代价 `1 / weight` 的最短路，用 Dijkstra 计算。

只在最大连通分量上选点。`k` 大于等于该分量的节点数时，返回整个分量。

两步：

1. **Farthest-first 初始化。** 用 `seed` 随机选第一个中心，之后每次选当前离已有中心最远的节点，直到有 `k` 个中心。
2. **受限 PAM。** 每一轮把节点分到最近的中心。对每个簇，只在簇内 strength 最高的 `kmedoids_m` 个节点里找替换中心。若替换能降低总距离，就换上。最多做 `kmedoids_max_iter` 轮，一轮没有任何替换就停止。

候选评估会缓存每个起点的最短路，避免对同一节点重复 Dijkstra。

### 1.4 参数

`gsemb_select_landmarks()` 与 `gsemb_fit()` 共用下面这些参数。`gsemb_fit()` 里方法名是 `landmark_method`，以免和嵌入方法 `method`（`svd` / `torch_autoencoder` / `set2gaussian_torch`）冲突。

| 参数 | 默认 | 谁用 | 含义 |
|---|---|---|---|
| `k` | 128 | 全部 | 要选的 landmark 个数，不会超过节点数。K-Medoids 在最大连通分量内再封顶。 |
| `method` / `landmark_method` | `"degree"` | 全部 | `"degree"`、`"random"`、`"betweenness"`、`"kmedoids"`。 |
| `seed` | 1 | `random`、`kmedoids` | 随机抽样或 farthest-first 第一个中心的随机种子。 |
| `betweenness_cutoff` | -1 | `betweenness` | 传给 `igraph::betweenness` 的路径长度上限，单位是边代价之和。`-1` 为完整介数。正数会丢掉长于该代价的路径。其他方法忽略。 |
| `kmedoids_m` | 50 | `kmedoids` | 每个簇里参与 PAM 替换的候选上限，按 strength 取前 M 个。越大越接近完整 PAM，也越慢。 |
| `kmedoids_max_iter` | 10 | `kmedoids` | PAM 最多轮数。提前没有改进会停。 |

`random` 不看权重，在全部节点里均匀抽 `k` 个。

拟合结果里的 `landmarks` 是选中的节点 ID，`landmark_method` 记录实际用的方法。

```r
fit <- gsemb_fit(  # edges 已按下文第 3 节清理
  ppi = edges,
  gene_sets = gene_sets,
  node1 = "protein1",
  node2 = "protein2",
  weight = "combined_score",
  landmark_method = "kmedoids",  # 或 "degree" / "betweenness" / "random"
  k = 128,
  seed = 1,
  kmedoids_m = 50,
  kmedoids_max_iter = 10
)
```

只选点、不做嵌入时：

```r
gsemb_select_landmarks(adj, k = 128, method = "betweenness")
```

---

## 2. 本次改动

老师确认保留现有加权实现后，接口收成四种方法，并接到主流程。这是一次**会改变旧结果的接口迁移**：旧主分支的 `method="degree"` 按邻居数量排序，新版同名方法按边权求和；在旧分析的记录中只有方法名、没有源码版本和图版本时，不能把两者当作同一种实验。

- 删掉无权对照：不再按邻居个数排 `degree`，不再用无权最短路介数，不再用 hop 距离做 K-Medoids。复现旧结果需使用旧代码及原图；使用新代码应重跑并记录图、方法、参数和种子。
- 去掉 `weighted_` 前缀。现在的 `degree` 就是原来的节点强度，`betweenness` 和 `kmedoids` 就是原来边代价为 `1/weight` 的版本。
- `gsemb_fit()`、`gsemb_compute_node_landmark_features()`、`gsemb_fit_set2gaussian_torch()` 增加 `landmark_method`，以及 `betweenness_cutoff`、`kmedoids_m`、`kmedoids_max_iter`。默认方法仍是 `degree`。
- `betweenness_cutoff` 默认从 `4` 改为 `-1`。在 STRING 高置信图上，`4` 大于所有路径代价，本来就不会截断；默认改成完整介数后，含义和结果一致。调用方仍可传入正数自行截断，参数名是 `betweenness_cutoff`，不是 igraph 的 `cutoff`。
- K-Medoids 的 C++ 只保留 Dijkstra 与按 strength 选候选，并缓存最短路。
- 上述包接口已进入 `main`（提交 `6eae7c8`）。研究脚本仍在分支 `feat/weighted-degree-landmarks`，没有并进 `main`。

### 2.1 实测耗时

图：STRING v12 人类 `9606.protein.links`，去掉自环和反向重复边后，`combined_score > 700`。16185 个节点，236000 条无向边。`k = 128`，`seed = 20260929`，`kmedoids_m = 50`，`kmedoids_max_iter = 10`。计时时介数 `cutoff = 4`；该图路径代价远小于 4，与现在的默认 `-1` 等价。时间为单次 `system.time` 的 elapsed，机器未记在表里。表中的旧方法名就是现在的 `degree` / `betweenness` / `kmedoids`。

来源：[原始计时表](benchmarks/weighted_landmarks_score700_timing_rectify0929.csv)（`-O2`，最短路按起点缓存）。

| 现在的方法 | 当时表中的名字 | elapsed（秒） |
|---|---|---|
| `degree` | `weighted_degree` | 0.003 |
| `betweenness` | `weighted_betweenness` | 179.739 |
| `kmedoids` | `weighted_kmedoids` | 39.622 |

同一张清理后的图、同一参数，在 `-O0` 且 K-Medoids 尚未缓存最短路时：`degree` 0.004 秒，`betweenness` 160.920 秒，`kmedoids` 1147.563 秒（[对照计时表](benchmarks/weighted_landmarks_score700_timing_rectify0929_O0.csv)）。两个运行同时改变了编译优化和最短路缓存，不能把约 19 分钟到约 40 秒的全部差异单独归因于缓存；要拆分贡献还需同编译配置下的消融计时。

更早一次计时用的是权重被加了两次的旧 high 图（边权约 1402–1998），且没有最短路缓存：`degree` 0.002 秒，`betweenness` 175.241 秒，`kmedoids` 1943.330 秒。那次不能和上表直接比选点名单。

这些数字只说明选点本身的耗时，不包含后面的 RWR、SVD 和富集。

### 2.2 三种方法的内存对照

本轮在 Linux、R 4.6.1 下用独立 R 进程和 `/usr/bin/time` 的最大常驻内存（peak RSS）测量。输入是本地历史 `test_results/high_conf_graph.rds`：16185 节点、236000 无向边，SHA-256 为 `ad5e60a070c4adbb65c4cb48559252c992dee9682cb77598ee09561f057d5a2d`。该缓存的反向边权曾相加到 1402–1998；内存比较使用相同的图规模和拓扑，**不能代替在新清理图上验证选点名单**。测试运行通过 `pkgload::load_all()` 装载包并使用调试编译，`k=128`、`seed=20260929`、`kmedoids_m=50`、`kmedoids_max_iter=10`、`betweenness_cutoff=-1`。原始读数见 [峰值内存表](benchmarks/landmark_peak_memory_20260929.csv)。

| 进程执行内容 | 峰值 RSS | 相比只加载图多用 | 说明 |
|---|---:|---:|---|
| 只加载包和图 | 约 219 MiB | 基线 | 包含 R 进程、已装载包和缓存图。 |
| `degree` | 约 219 MiB | 约 0.3 MiB | 只计算邻接矩阵行和，额外内存很小。 |
| `betweenness` | 约 274 MiB | 约 55 MiB | 构造 igraph 图并计算完整加权介数。 |
| `kmedoids` | 约 1382 MiB | 约 1163 MiB | 缓存各候选起点到所有节点的距离，随候选数增长。 |

表中是**进程峰值**，不是函数对象的单独大小；R 版本、系统、图规模、种子和候选数都会改变数值。这里 `kmedoids` 的调试构建进程耗时约 291 秒，不能直接与上面的 `-O2` 39.622 秒对比。当前最短路缓存没有内存上限：16185 节点时，每缓存一个起点约需 126 KiB 的距离数值；若缓存所有起点，距离数组本身约需 2.0 GiB。需要内存预算时，应在性能和峰值内存间做对照，再确定缓存上限。

---

## 3. `gsemb_clean_edges()`

`gsemb_build_graph()` 和 `gsemb_fit()` 只负责把边表变成邻接矩阵。无向时会做 `adj + t(adj)`，把每条边写成对称矩阵。如果边表里同一对蛋白两个方向都在，权重会被加两次。它们也不删自环，也不按分数过滤。

`gsemb_clean_edges()` 放在 `R/graph.R`，在构图之前可选地清洗边表。已经干净的自建网络可以不用它。

端点为 `NA` 或空字符串的行会先被跳过，并给出跳过的边数警告。介数和 K-Medoids 在选点入口检查非零边权是否为有限正数；零表示没有边，不要求输入权重属于 0–1 或 0–1000 的固定范围。

| 参数 | 默认 | 含义 |
|---|---|---|
| `node1`, `node2` | `"node1"`, `"node2"` | 端点列名。STRING 用 `"protein1"`, `"protein2"`。 |
| `weight` | `NULL` | 权重列。设置 `score_cutoff` 时必须给出。 |
| `score_cutoff` | `NULL` | 丢掉 `weight <= score_cutoff` 的边。`NULL` 表示不按分数过滤。STRING 高置信常用 `700`。 |
| `drop_self_loops` | `TRUE` | 丢掉两端点相同的行。 |
| `undirected_unique` | `TRUE` | 每个无向对只留一行，并把端点排成 `node1 < node2`。同一对有多行时，有权重则留最大权重，否则留第一行。 |

双向分数不一致时，这个函数只保留较大的那一行，不会把两边的权重相加。若实验需要合并权重，应在调用前自己处理。

```r
library(geneSetEmbedding)

raw <- read.table(
  "9606.protein.links.v12.0.txt",
  header = TRUE,
  sep = " ",
  stringsAsFactors = FALSE
)

edges <- gsemb_clean_edges(
  raw,
  node1 = "protein1",
  node2 = "protein2",
  weight = "combined_score",
  score_cutoff = 700
)

fit <- gsemb_fit(
  ppi = edges,
  gene_sets = gene_sets,
  node1 = "protein1",
  node2 = "protein2",
  weight = "combined_score",
  landmark_method = "degree"
)
```

自建、已经无自环且每条无向边只出现一次的网络，直接 `gsemb_fit()` 或 `gsemb_build_graph()` 即可。
