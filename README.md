# Clustering repeatedly observed event sequences by learning multi-marked point processes

## Referrence and Citation

Please cite following papers for utilizing the codes.

Lihao Yin, Ganggang Xu, Huiyan Sang, Yongtao Guan. [Row-clustering of a Point Process-valued Matrix](https://arxiv.org/pdf/2110.01207.pdf), NuerIPS 2021


## Introudction
Structured point process data harvested from various platforms poses new chal- lenges to the machine learning community. By imposing a matrix structure to repeatedly observed marked point processes, we propose a novel mixture model of multi-level marked point processes for identifying potential heterogeneity in the ob- served data. Specifically, we study a matrix whose entries are marked log-Gaussian Cox processes and cluster rows of such a matrix. An efficient semi-parametric Expectation-Solution (ES) algorithm combined with functional principal compo- nent analysis (FPCA) of point processes is proposed for model estimation. The effectiveness of the proposed framework is demonstrated through simulation studies and a real data analysis.


## Typical applications

## Code Implementation

### Learning the clusters

```
$ python3 main.py  --dataset_path   dataset_path   \
                   --label_path     label_path     \ 
                   --model          model          \
                   --nacounts       naacounts      \
                   --mdays          mdays          \
                   --event_types    event_types    \
                   --time_slot      time_slot      \
                   --nclusters      nclusters      \
                   --bwd            bwd            \
```
