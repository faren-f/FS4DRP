import decoupler as dc
import os
os.chdir("../Data/Read-TF_dR")
import pandas as pd

dorothea = dc.get_dorothea(organism='human', levels=['A','B','C'])


path_GE = '../Data/expresion_matrix.csv'

GE = pd.read_csv(path_GE, sep = ',', header = 0, index_col=0)
GE.head()
GE.shape

GE_norm = (GE - GE.mean())/GE.std()

tf_acts_gsea = dc.run_gsea(mat = GE_norm, net = dorothea, source = 'source', 
                           target = 'target', min_n = 5)

tf_acts_gsea = tf_acts_gsea[1]
