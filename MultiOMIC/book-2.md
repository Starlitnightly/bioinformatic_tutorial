## 2. GLUE多组学整合

文字版教程解释参照GLUE官方文档:https://scglue.readthedocs.io/zh_CN/latest/tutorials.html

但具体的代码参照我的文件，我对教程以及我们的数据进行了配套处理，分别是`glue-0.ipynb`,`glue-1.ipynb`,`glue-2.ipynb`三个文件。

但值得注意的是，我们这里得到的细胞只是整合后，还没进行配对，就是使得一个细胞同时具有两个组学层，配对的方法如下，我们需要使用具有X_glue层的rna与atac文件

```python
#读取数据
rna=sc.read('../cellanno/rna_anno.h5ad')
atac=sc.read('../cellanno/atac_anno.h5ad')

#提取GLUE层结果
rna_loc=pd.DataFrame(rna.obsm['X_glue'], index=rna.obs.index)
atac_loc=pd.DataFrame(atac.obsm['X_glue'], index=atac.obs.index)

#对GLUE层进行Pearson系数分析
import numpy as np
import gc
len1=(len(rna_loc)//5000)+1
p_pd=pd.DataFrame(columns=['rank_'+str(i) for i in range(50)])
n_pd=pd.DataFrame(columns=['rank_'+str(i) for i in range(50)])
for j in range(len1):
    c=pd.DataFrame()
    for i in range(len1):
        t1=rna_loc.iloc[5000*(i):5000*(i+1)]
        t2=atac_loc.iloc[5000*(j):5000*(j+1)]
        a=np.corrcoef(t1,t2)[len(t1):,0:len(t1)]
        b=pd.DataFrame(a,index=t2.index,columns=t1.index)  

        c=pd.concat([c,b],axis=1)
        del t1
        del t2
        del a
        del b
        gc.collect()
    for i in range(len(c)):
        t_c=c.iloc[i]
        p_pd.loc[t_c.name]=c.iloc[i].sort_values(ascending=False)[:50].values
        n_pd.loc[t_c.name]=c.iloc[i].sort_values(ascending=False)[:50].index.tolist()
    print('Now epoch is {}, {}/{}'.format(j,j*5000+len(c),len(atac_loc))) 
    del c
    gc.collect()
    
#寻找最近的细胞，其中depth的灵活调整可以使得配对成功的细胞数变大，同时精度有所下降
def find_neighbor_cell(p_pd,n_pd,depth=10):
    rubish_c=[]
    finish_c=[]
    for d in range(depth):
        p_pd=p_pd.loc[p_pd['rank_{}'.format(d)]>0.9]
        p_pd=p_pd.sort_values('rank_{}'.format(d),ascending=False)
        for i in p_pd.index:
            name=n_pd.loc[i,'rank_{}'.format(d)]
            if name not in rubish_c:
                finish_c.append(i)
                rubish_c.append(name)
            else:
                continue
        p_pd=p_pd.loc[~p_pd.index.isin(finish_c)]
        n_pd=n_pd.loc[~n_pd.index.isin(finish_c)]
    result=pd.DataFrame()
    result['omic_1']=rubish_c
    result['omic_2']=finish_c
    result.index=['cell_{}'.format(i) for i in range(len(result))]
    return result
  
res_pair=find_neighbor_cell(p_pd,n_pd,depth=20)
res_pair.head()
```

| omic_1 |             omic_2 |                        |
| -----: | -----------------: | ---------------------- |
| cell_0 |   AACTCAGCATGATCCA | neg-GCACGGTGTGCAAGAC-1 |
| cell_1 |   AAACGGGTCGGCGCAT | neg-TCTATTGAGAGGCAGG-1 |
| cell_2 |   GACGCGTAGACAAGCC | neg-GTGTCCTGTCATTGCA-1 |
| cell_3 | AAACCTGTCCCAACGG-1 | neg-AGTGTACCATGTGGGA-1 |
| cell_4 |   CATATTCCAAACCCAT | neg-TCAGCTCCAACGGGTA-1 |

## 3. 整体细胞分析

文字版教程解释参照Pyomic官方文档：https://pyomic.readthedocs.io/en/latest/Tutorials/t_cellanno.html

但具体的代码依旧参照我的文件，我对相关配置进行了修改处理，分别是`cellanno-1.ipynb`,`cellanno-2.ipynb`



