# 单细胞大模型专业中英文词汇表

这份词汇表面向单细胞基础模型、单细胞转录组分析和模型工程化落地。每个词条包含：

```text
英文 / 缩写
英文全称
中文
概念解释
实际意义与经验判断
```

使用建议：

```text
先掌握第 1-4 部分，能读懂单细胞数据和模型输入；
再掌握第 5-8 部分，能读懂论文和复现实验；
最后掌握第 9-12 部分，能做科研判断和工程落地。
```

## 1. 单细胞数据基础

| 英文/缩写 | 英文全称 | 中文 | 概念 | 实际意义与经验判断 |
|---|---|---|---|---|
| scRNA-seq | single-cell RNA sequencing | 单细胞 RNA 测序 | 在单细胞分辨率测量转录组表达。 | 所有单细胞大模型最常见输入。要先分清 raw counts、normalized counts、log-normalized data。 |
| snRNA-seq | single-nucleus RNA sequencing | 单核 RNA 测序 | 测量细胞核 RNA，常用于难以完整分离细胞的组织。 | 与 scRNA-seq 分布不同，核内转录本、低表达和细胞类型比例可能不同，模型迁移时要谨慎。 |
| UMI | unique molecular identifier | 唯一分子标识符 | 用随机 barcode 标记原始 RNA 分子以减少 PCR bias。 | UMI count 更接近分子数；非 UMI reads 与 UMI 数据不能简单混用。 |
| Count matrix | count matrix | 计数矩阵 | 细胞 × 基因的原始计数矩阵。 | 大模型 tokenizer 常要求 raw counts；把 log-normalized 数据误当 raw counts 是高频错误。 |
| Expression matrix | expression matrix | 表达矩阵 | 泛指基因表达矩阵，可为 raw、normalized 或 log-transformed。 | 写报告时不能只说 expression matrix，必须说明具体层。 |
| Cell | cell | 细胞 | 表达矩阵中的一个观测单位。 | 单细胞数据里的“cell”可能是细胞、细胞核或 barcode 液滴。 |
| Gene | gene | 基因 | 表达矩阵中的特征。 | 不同模型使用 gene symbol、Ensembl ID 或固定 vocabulary，映射失败会直接影响模型。 |
| Feature | feature | 特征 | 数据矩阵中的变量，可为基因、peak、蛋白、guide RNA 等。 | 多组学中 feature 不一定是 gene。 |
| Barcode | cell barcode | 细胞条形码 | 标记单个细胞/液滴的序列。 | barcode 不等于真实细胞；空液滴和 doublet 需要过滤。 |
| Library size | library size | 文库大小/总计数 | 每个细胞测到的总 counts。 | 影响表达强度；normalize_total 和 rank encoding 都在处理它带来的差异。 |
| Sparsity | sparsity | 稀疏性 | 表达矩阵中大量值为 0。 | 0 可能是真不表达，也可能是 dropout；模型和指标解释都要考虑。 |
| Dropout | dropout | 掉落/漏检 | 真实表达但未被测到。 | 单细胞数据噪声来源之一；不能把 0 值完全等同于无表达。 |
| Doublet | doublet | 双细胞/混合液滴 | 一个 barcode 捕获两个或多个细胞。 | 会形成假 cluster 或混合 marker，模型可能误学为新细胞类型。 |
| Ambient RNA | ambient RNA | 环境 RNA | 细胞破裂或背景导致的游离 RNA 污染。 | 常让不该表达的 marker 弱阳性；注释和模型输入都受影响。 |
| Batch | batch | 批次 | 技术批次、测序批次、样本处理批次等。 | batch 可能与疾病/样本强绑定，是数据泄漏和伪生物信号的核心来源。 |
| Donor | donor | 供体/个体 | 样本来源个体。 | 泛化评估应优先 donor-held-out，而不是随机细胞划分。 |
| Sample | sample | 样本 | 实验样本或生物样本。 | sample 与 batch、condition、donor 常混在一起，必须审计。 |
| Condition | condition | 条件 | 实验或生物条件，如 healthy/disease、treated/control。 | 不能让 condition 与 batch 完全混杂，否则模型可能只学 batch。 |
| Cell type | cell type | 细胞类型 | 稳定的细胞身份类别。 | 标签层级要统一；T cell、CD4 T、naive CD4 T 不是同一层级。 |
| Cell state | cell state | 细胞状态 | 同一细胞类型内的状态差异，如 activated、cycling、disease-associated。 | 大模型常更擅长连续状态表示，但标签和验证更难。 |
| Cell lineage | cell lineage | 细胞谱系 | 细胞发育和分化关系。 | 解释 embedding 或 trajectory 时要结合发育生物学证据。 |
| Marker gene | marker gene | 标志基因 | 某细胞类型/状态中相对特异表达的基因。 | 一个 marker 不能决定细胞类型；至少看 marker panel 和组织背景。 |
| Housekeeping gene | housekeeping gene | 管家基因 | 广泛表达、维持基础功能的基因。 | 高表达但不一定区分细胞状态；rank-value encoding 常试图降低其支配性。 |

## 2. 数据结构与文件格式

| 英文/缩写 | 英文全称 | 中文 | 概念 | 实际意义与经验判断 |
|---|---|---|---|---|
| AnnData | annotated data | 注释数据对象 | Python 单细胞分析常用对象。 | 读 h5ad 先查 `X/obs/var/obsm/layers/raw`，不要假设 `X` 是 raw counts。 |
| h5ad | HDF5 AnnData file | AnnData 文件 | AnnData 的磁盘格式。 | 单细胞大模型 demo 常以 h5ad 为输入；schema 不统一是落地难点。 |
| adata.X | AnnData X matrix | 主表达矩阵 | AnnData 中默认矩阵。 | 可能是 raw、log-normalized 或 scaled；使用前必须确认。 |
| adata.obs | observation metadata | 细胞元数据 | 每个细胞的标签、batch、donor 等信息。 | 数据泄漏检查主要看 obs。 |
| adata.var | variable metadata | 基因元数据 | 每个基因的 ID、symbol、HVG 标记等。 | Geneformer 需要 Ensembl ID；scGPT 需要 gene_col 匹配 vocabulary。 |
| adata.obsm | multidimensional observation annotation | 细胞多维表示 | PCA、UMAP、模型 embedding 常放这里。 | 建议命名为 `X_pca`、`X_scgpt`、`X_geneformer` 等。 |
| adata.layers | layers | 多层表达矩阵 | 保存 raw counts、normalized counts 等不同矩阵。 | 工程项目建议明确 `layers["counts"]`，避免覆盖原始矩阵。 |
| adata.raw | raw snapshot | 原始快照 | AnnData 中保存早期表达矩阵与 var。 | 不一定是真 raw counts；只是某一时刻快照。 |
| Loom | loom file | Loom 文件 | 单细胞矩阵文件格式之一。 | Geneformer tokenizer 支持 loom/h5ad；字段要求严格。 |
| MTX | matrix market format | 矩阵市场格式 | 10x 常见 sparse matrix 输出格式。 | 初始数据常为 `matrix.mtx + barcodes.tsv + features.tsv`。 |
| Sparse matrix | sparse matrix | 稀疏矩阵 | 只存非零值的矩阵格式。 | 大数据必须用稀疏矩阵；转 dense 容易爆内存。 |
| Gene symbol | gene symbol | 基因符号 | 如 IL7R、MS4A1。 | 直观但别名多、版本不稳定。 |
| Ensembl ID | Ensembl gene identifier | Ensembl 基因 ID | 如 ENSG00000168685。 | 更稳定；Geneformer 等模型强依赖。 |
| Vocabulary / Vocab | vocabulary | 词表 | 模型能识别的基因或 token 集合。 | gene coverage 是模型可用性的基础指标。 |
| Checkpoint | model checkpoint | 模型权重检查点 | 保存训练后参数的文件。 | 必须记录 checkpoint 名称、版本、来源、下载日期。 |
| Config / args.json | configuration file | 配置文件 | 保存模型结构和超参数。 | scGPT 等模型常需要 `args.json + best_model.pt + vocab.json` 同时匹配。 |

## 3. 质量控制、预处理与传统分析

| 英文/缩写 | 英文全称 | 中文 | 概念 | 实际意义与经验判断 |
|---|---|---|---|---|
| QC | quality control | 质量控制 | 过滤低质量细胞、低表达基因和异常样本。 | QC 阈值要看分布，不能机械照搬论文。 |
| n_genes_by_counts | number of genes by counts | 每细胞检测基因数 | 每个细胞非零基因数。 | 太低可能低质量，太高可能 doublet。 |
| total_counts | total counts | 总计数 | 每个细胞总表达量。 | 与 library size 类似，是测序深度 proxy。 |
| pct_counts_mt | percentage mitochondrial counts | 线粒体比例 | 线粒体基因 counts 占比。 | 高值可能表示细胞破损或应激，但组织差异很大。 |
| Ribosomal genes | ribosomal genes | 核糖体基因 | RPL/RPS 等。 | 高比例可能反映技术或生物状态，不能简单删除。 |
| Cell cycle score | cell cycle score | 细胞周期评分 | S/G2M 等细胞周期状态评分。 | 细胞周期可能是生物信号，也可能干扰聚类。 |
| Normalize total | normalize total counts | 总量归一化 | 把每个细胞总 counts 缩放到固定值。 | Scanpy 常用 `target_sum=1e4`；大模型 tokenizer 未必接受 normalized 数据。 |
| Log1p | log one plus transform | log(1+x) 转换 | 压缩高表达值动态范围。 | 常用于 PCA/传统分析；不要把 log 数据输入要求 raw counts 的 tokenizer。 |
| Scaling | feature scaling | 标准化/缩放 | 让基因均值为 0、方差为 1。 | PCA/线性模型常用；稀疏矩阵 dense 化要小心。 |
| HVG | highly variable genes | 高变基因 | 选择跨细胞变化大的基因。 | 传统分析常用；Geneformer 明确不应先做 HVG 后 tokenization。 |
| PCA | principal component analysis | 主成分分析 | 线性降维。 | 必须保留作为大模型 embedding 的最低 baseline。 |
| kNN graph | k-nearest neighbor graph | k 近邻图 | 用细胞相似性构建邻接图。 | Leiden/UMAP 常基于它；n_neighbors 会影响结果。 |
| UMAP | uniform manifold approximation and projection | UMAP 降维图 | 非线性二维可视化方法。 | 只能辅助观察，不能作为唯一证据。 |
| t-SNE | t-distributed stochastic neighbor embedding | t-SNE 降维图 | 非线性可视化方法。 | 局部结构强，整体距离更不能过度解释。 |
| Leiden | Leiden clustering | Leiden 聚类 | 图社区发现算法。 | resolution 影响 cluster 数；不能把 cluster 直接等同真实细胞类型。 |
| Louvain | Louvain clustering | Louvain 聚类 | 早期常用图聚类。 | 现多用 Leiden，但旧论文常见。 |
| DE genes | differentially expressed genes | 差异表达基因 | 条件或群体间表达显著差异的基因。 | 扰动预测常看 DE overlap；统计显著不等于生物重要。 |
| Pseudobulk | pseudobulk | 伪 bulk | 按样本/条件聚合单细胞 counts 后做 bulk-like 分析。 | 疾病差异分析常比逐细胞统计更稳健。 |
| Ambient correction | ambient RNA correction | 环境 RNA 校正 | 去除背景 RNA 污染。 | 对低表达 marker 和罕见细胞很关键。 |
| Doublet detection | doublet detection | 双细胞检测 | 识别混合 barcode。 | doublet 会误导模型学习假细胞状态。 |

## 4. 生物统计与评价指标

| 英文/缩写 | 英文全称 | 中文 | 概念 | 实际意义与经验判断 |
|---|---|---|---|---|
| Accuracy | accuracy | 准确率 | 总体预测正确比例。 | 类别不平衡时容易虚高。 |
| Precision | precision | 精确率 | 预测为某类中真正为该类的比例。 | 适合看误报。 |
| Recall | recall | 召回率 | 真实某类中被找出的比例。 | 稀有细胞类型尤其要看 recall。 |
| F1 score | F1 score | F1 分数 | precision 与 recall 的调和平均。 | 比 accuracy 更适合类别不平衡任务。 |
| Macro-F1 | macro-averaged F1 | 宏平均 F1 | 对每类 F1 等权平均。 | 细胞注释必须报告；能暴露小类失败。 |
| Micro-F1 | micro-averaged F1 | 微平均 F1 | 按样本汇总后计算 F1。 | 多数类占主导，接近 accuracy。 |
| Confusion matrix | confusion matrix | 混淆矩阵 | 真实标签与预测标签交叉表。 | 失败分析核心图。 |
| ARI | adjusted rand index | 调整兰德指数 | 聚类与真实标签一致性指标。 | 聚类评估常用；对随机一致性做校正。 |
| NMI | normalized mutual information | 归一化互信息 | 聚类与标签信息共享程度。 | 与 ARI 搭配看。 |
| Silhouette score | silhouette score | 轮廓系数 | 衡量同类近、异类远。 | `silhouette(cell type)` 和 `silhouette(batch)` 含义相反，解释要小心。 |
| kBET | k-nearest neighbor batch effect test | kNN 批次效应检验 | 检查局部邻域 batch 混合是否符合全局比例。 | 批次整合常用，但对参数和样本比例敏感。 |
| LISI | local inverse Simpson's index | 局部逆 Simpson 指数 | 局部邻域标签多样性指标。 | iLISI 看 batch mixing，cLISI 看 cell type purity。 |
| iLISI | integration LISI | 整合 LISI | 衡量 batch 混合。 | 越高常表示 batch 混合越好，但不要混掉生物状态。 |
| cLISI | cell-type LISI | 细胞类型 LISI | 衡量 cell type 是否混杂。 | 细胞类型应保持区分，cLISI 过高可能过整合。 |
| Pearson correlation | Pearson correlation | 皮尔逊相关 | 线性相关。 | 扰动表达预测常用；可能掩盖关键基因方向错误。 |
| Spearman correlation | Spearman rank correlation | 斯皮尔曼相关 | 排名相关。 | 对非线性和极值更稳健。 |
| AUROC | area under ROC curve | ROC 曲线下面积 | 二分类排序能力。 | 类别不平衡时还需看 AUPRC。 |
| AUPRC | area under precision-recall curve | PR 曲线下面积 | 不平衡二分类评估。 | 稀有事件预测更有意义。 |
| Baseline | baseline model | 基线模型 | 用于比较的简单或经典方法。 | 大模型必须超过合理 baseline 才能说有额外价值。 |
| Benchmark | benchmark | 基准评测 | 固定数据、任务、指标的系统比较。 | 公平 benchmark 比单个高分更重要。 |
| Held-out split | held-out split | 留出划分 | 训练时不使用的测试划分。 | donor-held-out 比 random cell split 更可信。 |
| Cross-validation | cross-validation | 交叉验证 | 多次划分训练/测试评估稳定性。 | 小数据项目建议使用，但注意按 donor 分组。 |
| Negative control | negative control | 阴性对照 | 理论上不应产生正结果的对照。 | 判断模型是否走捷径。 |
| Sanity check | sanity check | 合理性检查 | 检查结果是否符合基本预期。 | 如打乱标签后性能应下降。 |

## 5. 深度学习与 Transformer

| 英文/缩写 | 英文全称 | 中文 | 概念 | 实际意义与经验判断 |
|---|---|---|---|---|
| DL | deep learning | 深度学习 | 用多层神经网络学习表示。 | 单细胞大模型是深度学习在转录组上的扩展。 |
| Neural network | neural network | 神经网络 | 由可学习参数组成的函数。 | 不要把复杂网络自动等同可解释机制。 |
| Parameter | model parameter | 模型参数 | 训练中更新的权重。 | 参数量越大不等于生物学能力越强。 |
| Hyperparameter | hyperparameter | 超参数 | 人为设置的训练/模型参数。 | 学习率、batch size、input length 会显著影响结果。 |
| Loss | loss function | 损失函数 | 衡量预测与目标差距。 | 训练目标决定模型学什么。 |
| Epoch | epoch | 训练轮次 | 遍历训练数据一遍。 | fine-tuning 过多可能过拟合。 |
| Batch size | batch size | 批大小 | 一次训练/推理的样本数。 | 显存不足先减 batch size。 |
| Embedding | embedding | 嵌入向量 | 把离散对象变成连续向量。 | gene/cell embedding 是下游分析的核心中间产物。 |
| Token | token | 词元 | 模型输入的最小单位。 | 单细胞中 token 可为基因、表达值、metadata、扰动条件。 |
| Hidden state | hidden state | 隐状态 | 模型中间层输出表示。 | 不同层 embedding 含义不同；最后层不一定最适合泛化。 |
| Transformer | transformer | Transformer 模型 | 基于 self-attention 的深度模型架构。 | 单细胞模型常用，但基因没有天然语言顺序。 |
| Attention | attention | 注意力机制 | 让 token 根据相关性聚合其他 token 信息。 | attention 权重不等于因果调控。 |
| Self-attention | self-attention | 自注意力 | 同一序列内部 token 相互关注。 | 用于学习基因共表达上下文。 |
| Multi-head attention | multi-head attention | 多头注意力 | 多组 attention 并行学习不同关系。 | 不同 head 不一定都有生物意义。 |
| Encoder | encoder | 编码器 | 把输入序列编码成表示。 | scBERT、Geneformer 多为 encoder 路线。 |
| Decoder | decoder | 解码器 | 根据上下文生成输出。 | 生成式任务常用。 |
| Encoder-decoder | encoder-decoder | 编码器-解码器 | 先编码再生成/重建。 | 多组学转换、表达重建常见。 |
| CLS token | classification token | 分类词元 | 用于汇聚整个序列信息的特殊 token。 | scGPT/Geneformer V2 中常用于 cell embedding。 |
| PAD token | padding token | 填充词元 | 填充序列到同一长度。 | mask/padding 错误会污染模型表示。 |
| Mask token | mask token | 遮蔽词元 | 训练中遮住部分输入让模型预测。 | masked modeling 的基础。 |
| MLM | masked language modeling | 遮蔽语言建模 | 遮住 token 并预测。 | scBERT/Geneformer 类模型常用思想。 |
| Autoregressive modeling | autoregressive modeling | 自回归建模 | 逐步预测下一个 token。 | 单细胞中基因无天然顺序，需明确序列顺序。 |
| Generative pretraining | generative pretraining | 生成式预训练 | 学习生成/补全数据模式的预训练。 | scGPT 核心概念，不等于生成假细胞就是最终目的。 |
| Fine-tuning | fine-tuning | 微调 | 在预训练模型上继续用任务数据训练。 | 使用目标标签更新参数，不是 zero-shot。 |
| Feature extraction | feature extraction | 特征抽取 | 不更新模型，只抽 embedding 做下游任务。 | 课程项目最稳的上手方式。 |
| Zero-shot | zero-shot learning | 零样本学习 | 不用目标任务标签直接应用模型。 | 必须检查预训练语料重叠，不能滥称泛化。 |
| Few-shot | few-shot learning | 少样本学习 | 少量标签适配模型。 | 小样本疾病任务常见，但 split 要严格。 |
| Transfer learning | transfer learning | 迁移学习 | 从大数据学表示迁移到小任务。 | Geneformer 的核心范式。 |
| Overfitting | overfitting | 过拟合 | 模型记住训练数据，泛化差。 | fine-tuning 小数据时常见。 |
| Domain shift | domain shift | 域偏移 | 训练与测试数据分布不同。 | 跨组织、跨平台、跨物种时核心挑战。 |
| OOD | out-of-distribution | 分布外 | 测试数据超出训练分布。 | unknown cell type、罕见疾病状态常属 OOD。 |

## 6. 单细胞大模型输入表示

| 英文/缩写 | 英文全称 | 中文 | 概念 | 实际意义与经验判断 |
|---|---|---|---|---|
| Tokenization | tokenization | 词元化 | 把数据转为 token 序列。 | 单细胞 tokenization 是建模假设，不是小预处理。 |
| Gene token | gene token | 基因词元 | 每个基因作为 token。 | 依赖 vocabulary 和 ID 映射。 |
| Value token | expression value token | 表达值词元 | 表达强度离散化后的 token。 | value binning 的基础。 |
| Value binning | value binning | 表达值分箱 | 把连续表达值分成有限档。 | 稳健但损失精细强度。 |
| Value projection | value projection | 表达值投影 | 把连续表达值投影成向量。 | 保留强度但对归一化敏感。 |
| Rank-value encoding | rank-value encoding | 排名-表达编码 | 根据细胞内基因表达优先级排序。 | Geneformer 核心；抗深度差异但损失绝对强度。 |
| Gene ranking | gene ranking | 基因排序 | 按表达或重要性对基因排序。 | top-k 太小会丢信号，太大可能引噪声。 |
| Top-k genes | top-k genes | 前 k 个基因 | 每细胞保留排名靠前基因。 | 影响模型输入长度和信息保留。 |
| Metadata token | metadata token | 元数据词元 | 组织、物种、batch、condition 等作为 token。 | 可能增强条件建模，也可能标签泄漏。 |
| Perturbation token | perturbation token | 扰动词元 | 药物、基因敲除等扰动条件 token。 | 扰动预测中常见。 |
| Species token | species token | 物种词元 | 表示 human/mouse 等物种。 | 跨物种模型需要，但 ortholog 映射仍关键。 |
| Tissue token | tissue token | 组织词元 | 表示组织来源。 | 可能帮助模型，也可能提示答案。 |
| Batch token | batch token | 批次词元 | 表示技术批次。 | 易造成泄漏；使用时必须审计。 |
| Gene embedding | gene embedding | 基因嵌入 | 基因的向量表示。 | 可做基因相似性假设，但不能直接当调控关系。 |
| Cell embedding | cell embedding | 细胞嵌入 | 细胞级向量表示。 | 聚类、注释、参考映射的核心。 |
| CLS embedding | CLS embedding | CLS 嵌入 | CLS token 的 hidden state。 | 常用于 cell-level 表示。 |
| Mean pooling | mean pooling | 平均池化 | 对 token hidden states 求平均。 | 简单稳健，但可能稀释关键基因。 |
| Attention pooling | attention pooling | 注意力池化 | 学习加权聚合 token。 | 更灵活，但解释更复杂。 |
| Gene coverage | gene coverage | 基因覆盖率 | 数据基因匹配模型词表比例。 | 工程落地必须报告；低覆盖会削弱结果可信度。 |

## 7. 代表模型与缩写

| 英文/缩写 | 英文全称 | 中文 | 概念 | 实际意义与经验判断 |
|---|---|---|---|---|
| scFM | single-cell foundation model | 单细胞基础模型 | 面向单细胞数据的大规模预训练模型。 | 不同模型任务边界差异很大，不能只看“foundation”。 |
| FM | foundation model | 基础模型 | 在大规模数据上预训练，可迁移到多任务的模型。 | 迁移能力必须用公平 benchmark 证明。 |
| scBERT | single-cell BERT | 单细胞 BERT 模型 | BERT-like encoder 用于细胞类型注释。 | 学习监督微调和 unknown detection 的好案例。 |
| Geneformer | Geneformer | Geneformer 模型 | rank-value encoding + masked gene prediction 的网络生物学模型。 | raw counts、Ensembl ID、n_counts 是关键输入要求。 |
| Genecorpus | Genecorpus | Geneformer 预训练语料 | Geneformer 的大规模单细胞语料。 | 预训练域决定模型泛化边界。 |
| scGPT | single-cell GPT | 单细胞 GPT 模型 | 生成式单细胞基础模型，支持多任务。 | 官方 tutorial 丰富，适合 embedding/reference mapping 入门。 |
| scFoundation | scFoundation | 单细胞基础模型 scFoundation | 大规模人类单细胞转录组基础模型。 | 适合人类细胞/基因表示，但需检查词表和任务域。 |
| UCE | universal cell embeddings | 通用细胞嵌入 | 强调 zero-shot universal cell embedding。 | 适合作为快速 embedding benchmark 对象。 |
| CellFM | Cell Foundation Model | 细胞基础模型 CellFM | 约 1 亿人类细胞预训练的大规模模型。 | 参数大、工程门槛高，需看成本收益。 |
| scLong | single-cell Long-context model | 长上下文单细胞模型 | 十亿参数级，强调长程基因上下文。 | 适合全基因长程依赖问题，但成本高。 |
| Performer | Performer | Performer 模型 | 高效 attention 变体。 | scBERT 使用，适合长序列降低 attention 成本。 |
| RetNet | retention network | RetNet/保留网络 | 一类高效序列建模架构。 | CellFM 相关架构关键词，理解为高效长序列替代路线。 |
| FAISS | Facebook AI Similarity Search | 向量相似搜索库 | 大规模最近邻搜索工具。 | scGPT reference mapping 等可用它加速。 |
| Model zoo | model zoo | 模型库 | 官方提供的多个 checkpoint 集合。 | checkpoint 选择会显著影响结果。 |

## 8. 单细胞大模型任务

| 英文/缩写 | 英文全称 | 中文 | 概念 | 实际意义与经验判断 |
|---|---|---|---|---|
| Cell type annotation | cell type annotation | 细胞类型注释 | 给细胞预测类型标签。 | 必须报告 macro-F1 和混淆矩阵。 |
| Reference mapping | reference mapping | 参考映射 | 把 query 映射到 reference embedding/atlas。 | 常用 embedding + nearest neighbor，比 fine-tuning 更轻。 |
| Batch integration | batch integration | 批次整合 | 减少技术批次差异。 | 不能只追求 batch 混合，还要保留细胞类型。 |
| Data integration | data integration | 数据整合 | 合并多个数据集/批次/平台。 | 标签、平台、组织混杂时最容易误判。 |
| Multi-omics integration | multi-omics integration | 多组学整合 | 整合 RNA、ATAC、蛋白等模态。 | 特征空间不同，不能简单拼接。 |
| Perturbation prediction | perturbation prediction | 扰动预测 | 预测基因/药物扰动后的表达变化。 | 不能等同因果证明。 |
| In silico perturbation | in silico perturbation | 计算机模拟扰动 | 在模型中模拟删除/上调/干预基因。 | 候选生成工具，需要实验验证。 |
| In silico treatment | in silico treatment | 计算治疗模拟 | 模拟让疾病状态向健康状态移动的扰动。 | 靶点优先级而非临床结论。 |
| GRN | gene regulatory network | 基因调控网络 | 基因间调控关系网络。 | attention/embedding 只能提供假设，不是直接 GRN。 |
| Trajectory inference | trajectory inference | 轨迹推断 | 推断细胞发育或状态连续变化路径。 | embedding 好坏会影响轨迹，需生物背景验证。 |
| Pseudotime | pseudotime | 伪时间 | 表示细胞沿轨迹的相对进程。 | 不是实际时间，不能过度解释。 |
| Cell state classification | cell state classification | 细胞状态分类 | 区分 healthy/disease、activated/resting 等状态。 | 必须按 donor 或 sample split 才可信。 |
| Disease classification | disease classification | 疾病分类 | 基于细胞或样本表达预测疾病状态。 | 细胞级 accuracy 常被样本泄漏夸大。 |
| Target discovery | therapeutic target discovery | 治疗靶点发现 | 寻找可能干预疾病的基因/通路。 | 模型只能提出候选，需湿实验验证。 |
| Expression imputation | expression imputation | 表达补全 | 预测缺失或 dropout 表达。 | 可能引入模型幻觉，不能替代真实测量。 |
| Denoising | denoising | 去噪 | 从噪声表达中恢复信号。 | 去噪可能抹掉真实稀有状态。 |
| Atlas building | atlas building | 图谱构建 | 构建大规模参考细胞图谱。 | 模型可帮助映射，但 atlas 标签质量决定上限。 |

## 9. 数据泄漏、复现与公平性

| 英文/缩写 | 英文全称 | 中文 | 概念 | 实际意义与经验判断 |
|---|---|---|---|---|
| Data leakage | data leakage | 数据泄漏 | 测试信息进入训练或调参。 | foundation model 项目最危险问题之一。 |
| Pretraining overlap | pretraining overlap | 预训练重叠 | 测试数据可能出现在预训练语料中。 | 不能排除时，不要声称严格 zero-shot。 |
| Label leakage | label leakage | 标签泄漏 | 输入中包含直接或间接标签信息。 | metadata token 和 batch/condition 绑定常见。 |
| Batch confounding | batch confounding | 批次混杂 | batch 与生物条件强相关。 | 模型可能学技术差异而非生物学。 |
| Donor leakage | donor leakage | 供体泄漏 | 同一 donor 同时出现在 train/test。 | 疾病分类和注释任务会虚高。 |
| Sample leakage | sample leakage | 样本泄漏 | 同一样本细胞跨 train/test。 | 细胞级随机 split 的常见问题。 |
| Train/test split | train/test split | 训练/测试划分 | 分离训练和测试数据。 | split 方式决定结论强度。 |
| Validation set | validation set | 验证集 | 用于调参和早停的数据。 | 不要用测试集调 threshold。 |
| Test set | test set | 测试集 | 只用于最终报告。 | 多次看测试集会形成隐性泄漏。 |
| Donor-held-out | donor-held-out | 供体留出 | 按 donor 留出测试。 | 比随机细胞划分更接近真实泛化。 |
| Batch-held-out | batch-held-out | 批次留出 | 按 batch 留出测试。 | 测试跨批次鲁棒性。 |
| Study-held-out | study-held-out | 研究留出 | 按数据集/研究留出。 | 最强泛化测试之一。 |
| Reproducibility | reproducibility | 可复现性 | 其他人能重跑并得到相近结果。 | 记录环境、seed、checkpoint、数据版本。 |
| Random seed | random seed | 随机种子 | 控制随机性的数值。 | 单一 seed 不够稳，关键结论建议多 seed。 |
| Ablation | ablation study | 消融实验 | 去掉某组件观察影响。 | 判断模型组件是否真的有贡献。 |
| Failure case | failure case | 失败案例 | 模型表现差的样本/类别。 | 比只展示成功图更有价值。 |
| Error analysis | error analysis | 错误分析 | 系统分析预测错误。 | 混淆矩阵、rare class、batch-specific error 必看。 |
| Model card | model card | 模型说明卡 | 记录模型训练、用途、限制。 | 工程落地应为每个 checkpoint 建 model card。 |
| Data card | data card | 数据说明卡 | 记录数据来源、处理、限制。 | 防止数据含义被误用。 |

## 10. 可解释性、因果与生物学判断

| 英文/缩写 | 英文全称 | 中文 | 概念 | 实际意义与经验判断 |
|---|---|---|---|---|
| Interpretability | interpretability | 可解释性 | 理解模型为何给出结果。 | 解释性不是自动等于生物机制。 |
| Explainability | explainability | 可解释/可说明性 | 对模型输出提供解释。 | 常与 interpretability 混用。 |
| Attribution | feature attribution | 特征归因 | 估计输入特征对输出贡献。 | 对输入尺度和模型不稳定性敏感。 |
| Integrated gradients | integrated gradients | 积分梯度 | 一种梯度归因方法。 | 可用于基因贡献，但需 baseline 选择。 |
| Occlusion | occlusion analysis | 遮挡分析 | 移除某特征看预测变化。 | 类似 in silico gene deletion。 |
| Counterfactual | counterfactual | 反事实 | 如果改变某因素，结果如何变化。 | 扰动预测的核心思路，但不等于真实因果。 |
| Causality | causality | 因果性 | 一个因素改变导致另一个因素改变。 | 单细胞观测数据主要是相关，不是因果。 |
| Correlation | correlation | 相关性 | 两变量共同变化。 | 共表达不等于调控。 |
| Co-expression | co-expression | 共表达 | 多基因表达模式相似或共同出现。 | 大模型常学习它，但需区别调控关系。 |
| Pathway enrichment | pathway enrichment | 通路富集 | 判断基因集是否富集某通路。 | 常用于解释，但容易后验讲故事。 |
| GO | gene ontology | 基因本体 | 基因功能、过程、组分注释体系。 | scLong 等可能使用 GO 图，跨物种需谨慎。 |
| TF | transcription factor | 转录因子 | 调控基因转录的蛋白。 | 低表达 TF 容易被过滤或 rank 表示低估/重排。 |
| GRN inference | gene regulatory network inference | 基因调控网络推断 | 从数据推断调控关系。 | 需要独立验证，不能只靠 attention。 |
| Marker validation | marker validation | 标志基因验证 | 用已知 marker 支持注释。 | 新状态必须结合 marker panel 和文献。 |
| External validation | external validation | 外部验证 | 独立数据或实验验证。 | 从模型发现到生物结论的关键。 |
| Wet-lab validation | wet-lab validation | 湿实验验证 | CRISPR、药物、流式、染色等实验验证。 | 靶点发现类项目最终必需。 |
| Biological plausibility | biological plausibility | 生物学合理性 | 是否符合已有知识。 | 只能支持假设，不是决定性证据。 |

## 11. 多组学、空间与扰动

| 英文/缩写 | 英文全称 | 中文 | 概念 | 实际意义与经验判断 |
|---|---|---|---|---|
| scATAC-seq | single-cell ATAC sequencing | 单细胞染色质开放性测序 | 测量染色质可及性。 | 与 RNA 特征空间不同，整合需特别方法。 |
| CITE-seq | cellular indexing of transcriptomes and epitopes by sequencing | CITE-seq | 同时测 RNA 和表面蛋白。 | 多组学对齐和标签验证很有用。 |
| Multiome | multiome | 多组学联合测量 | 同一细胞多模态数据。 | matched multiome 比 unmatched integration 更可靠。 |
| Spatial transcriptomics | spatial transcriptomics | 空间转录组 | 保留空间位置的表达测量。 | 单细胞大模型还需处理邻域和组织结构。 |
| Visium | Visium | 10x 空间转录组平台 | spot-level 空间表达。 | spot 可能含多个细胞，不等于单细胞。 |
| Xenium | Xenium | 高分辨原位平台 | 单细胞/亚细胞空间表达。 | 基因 panel 限制影响模型输入。 |
| Perturb-seq | Perturb-seq | 扰动单细胞测序 | CRISPR 扰动结合 scRNA-seq。 | 扰动预测的重要真实验证数据。 |
| CRISPR | clustered regularly interspaced short palindromic repeats | CRISPR | 基因编辑/扰动技术。 | 模型靶点假设常需 CRISPR 验证。 |
| KO | knockout | 基因敲除 | 使基因失活。 | 扰动类型之一。 |
| KD | knockdown | 基因敲低 | 降低基因表达。 | 与 KO 效果不同，不能混用。 |
| OE | overexpression | 过表达 | 提高基因表达。 | Geneformer/scGPT 可能模拟 overexpress。 |
| Guide RNA / gRNA | guide RNA | 向导 RNA | CRISPR 定位序列。 | Perturb-seq 中 guide assignment 质量影响标签。 |
| DE overlap | differential expression overlap | 差异基因重叠 | 预测与真实 DE genes 的重叠。 | 扰动预测比全基因相关更关注功能变化。 |
| Direction accuracy | direction accuracy | 方向准确率 | 表达变化方向是否预测正确。 | 靶点和扰动任务非常重要。 |
| Dose response | dose response | 剂量反应 | 不同剂量导致的响应差异。 | 大多数模型对剂量泛化仍有限。 |
| Time course | time course | 时间序列 | 不同时间点观测。 | 细胞状态动态建模难点。 |

## 12. 工程化、MLOps 与落地

| 英文/缩写 | 英文全称 | 中文 | 概念 | 实际意义与经验判断 |
|---|---|---|---|---|
| Pipeline | pipeline | 流水线 | 数据处理到模型输出的自动流程。 | 工程落地必须标准化。 |
| MLOps | machine learning operations | 机器学习运维 | 模型部署、监控、版本管理。 | 单细胞模型落地需要数据和模型双版本控制。 |
| Data validation | data validation | 数据校验 | 检查输入 schema、字段、矩阵类型。 | 防止 raw/log 层混用。 |
| Schema | schema | 数据模式 | 文件字段和格式规范。 | h5ad 工程化必须定义 obs/var/layers 必填项。 |
| Model registry | model registry | 模型注册表 | 管理模型版本和元数据。 | 记录 checkpoint、vocab、license、硬件需求。 |
| Versioning | versioning | 版本管理 | 管理数据、代码、模型版本。 | 不记录版本的模型结果不可复现。 |
| Git commit | git commit | Git 提交版本 | 代码版本标识。 | 报告中应记录官方仓库 commit 或 package version。 |
| Environment | environment | 运行环境 | Python、CUDA、PyTorch、包版本。 | flash-attn/CUDA 问题常导致 scGPT 等安装失败。 |
| CUDA | compute unified device architecture | CUDA | NVIDIA GPU 计算平台。 | 许多大模型依赖 GPU/CUDA 版本匹配。 |
| GPU memory | GPU memory | 显存 | GPU 可用内存。 | batch size、input length、模型规模受它限制。 |
| Runtime | runtime | 运行时间 | 模型运行耗时。 | 工程可用性指标，不是附属信息。 |
| Throughput | throughput | 吞吐量 | 单位时间处理细胞数。 | 大规模 atlas embedding 服务必须关注。 |
| Quantization | quantization | 量化 | 降低参数精度以减小模型和加速。 | 未来落地方向，但可能影响精度。 |
| Distillation | knowledge distillation | 知识蒸馏 | 用大模型指导小模型。 | 适合把大模型能力转为轻量部署。 |
| LoRA | low-rank adaptation | 低秩适配 | 参数高效微调方法。 | 单细胞模型低资源 fine-tuning 的潜在方向。 |
| PEFT | parameter-efficient fine-tuning | 参数高效微调 | 只训练少量参数。 | 降低显存和过拟合风险。 |
| Monitoring | monitoring | 监控 | 跟踪模型运行和数据漂移。 | 新数据 gene coverage、batch drift、失败率都要监控。 |
| Data drift | data drift | 数据漂移 | 新数据分布与历史不同。 | 医院/平台/组织变化会导致模型退化。 |
| Audit trail | audit trail | 审计轨迹 | 保留每步记录。 | 临床和科研复现都需要。 |
| Human-in-the-loop | human in the loop | 人在回路 | 模型输出由专家审查。 | 单细胞注释和生物发现不应完全自动化。 |
| License | license | 许可证 | 代码/模型/数据使用许可。 | 工业落地必须检查。 |
| Privacy | privacy | 隐私 | 样本和患者信息保护。 | 单细胞数据可能含遗传和疾病敏感信息。 |
| Federated learning | federated learning | 联邦学习 | 数据不出机构的联合训练。 | 医疗多中心单细胞建模的重要方向。 |

## 13. 老手判断清单

读任何单细胞大模型论文或项目，按下面清单问一遍：

```text
1. 输入是什么？raw counts、log expression、ranking、binning 还是 projection？
2. 基因 ID 是 gene symbol 还是 Ensembl ID？gene coverage 多少？
3. 模型训练数据覆盖我的物种、组织、疾病和平台吗？
4. 这是 zero-shot、feature extraction、reference mapping 还是 fine-tuning？
5. train/test split 是 random cell、donor-held-out、batch-held-out 还是 study-held-out？
6. 是否可能有预训练语料重叠？
7. 是否有合理 baseline？
8. 是否报告 macro-F1、混淆矩阵、batch mixing 和失败案例？
9. 是否只展示 UMAP？
10. 模型关注的是 marker、batch、cell cycle、mitochondrial genes 还是疾病通路？
11. attention/attribution 是否被过度解释为因果？
12. 扰动预测是否有外部实验或独立数据验证？
13. 运行时间、显存、失败细胞数是否可接受？
14. 结论是否明确写出不能证明什么？
```

经验上，能把这些问题答清楚，比会背模型名字重要得多。
