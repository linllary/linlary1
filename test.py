import pandas as pd


#First Part
A1=pd.read_excel('Python_testdata.xlsx',engine='openpyxl',sheet_name='test01_SV')
A2=pd.read_excel('Python_testdata.xlsx',engine='openpyxl',sheet_name='test02_SV')
B1=pd.read_excel('Python_testdata.xlsx',engine='openpyxl',sheet_name='test01_SNV')
B2=pd.read_excel('Python_testdata.xlsx',engine='openpyxl',sheet_name='test02_SNV')
B3=pd.read_excel('Python_testdata.xlsx',engine='openpyxl',sheet_name='output')

def counts(A1,A2,N):
 A1[A2.columns[N]]=A2[A2.columns[N]]
 A1['Counts']=2

 k=0
 for f1,f2 in zip(A1[A1.columns[N]],A1.columns[N-1]):
    #print(f1,f2)
    if f1=='.' or f2=='.':
       A1['Counts'].values[k]=1
    k+=1
 return A1

A1=counts(A1,A2,A1.shape[1]-1)
B1=counts(B1,B2,B1.shape[1]-1)
A1['Candidate_gene_filter']='.'
B1['Candidate_gene_filter']='.'

A1.rename({A1.columns[0]: B1.columns[0], A1.columns[1]: B1.columns[1], A1.columns[2]: B1.columns[2], A1.columns[13]: 'test01', A1.columns[14]: 'test02'}, axis=1, inplace=True)
A1.rename({'Gene name':'Gene_name','location1': 'Func.refGene/Location1','location2' : 'ExonicFunc.refGene/location2','GD_ID':'Variant_ID','GD_AF':'GD_AF_EAS','AnnotSV ranking':'ACMG_classes/AnnotSV_ranking' }, axis=1, inplace=True)
B1.rename({'Gene.refGene':'Gene_name', 'AAChange.refGene':'AAChange', 'avsnp150':'Variant_ID' , 'Func.refGene' : 'Func.refGene/Location1','ExonicFunc.refGene' : 'ExonicFunc.refGene/location2', 'AF_eas.1':'GD_AF_EAS',}, axis=1, inplace=True)

C=pd.concat([B1,A1], ignore_index=True)

C=C.sort_values(by=['Start'])

C.rename({'SV length' : 'SV_length', 'SV type':'SV_type'},axis=1,inplace=True )

k=0
for f in B3.columns:
    if k==0:
       D=C[f]
    else:
       D=pd.concat([D,C[f]],axis=1, ignore_index=True)
    k+=1

for f1,f2 in zip(D.columns,B3.columns):
    D.rename({f1:f2},axis=1,inplace=True )
k=0
for f1,f2,f3 in zip(D['test01'],D['test02'],D['Counts']):
    if f1=='.' or f2=='.':
       D['Counts'].values[k]=1
    else:
       D['Counts'].values[k]=2
    k+=1
k=0
for f1,f2 in zip(D['Gene_name'],D['Candidate_gene_filter']):
    if f1 =='DVL1' or f1=='SDF4':
       D['Candidate_gene_filter'].values[k]='+'
    k+=1
#print(D)
D.to_csv('output1.csv')


#second part


S1=pd.read_csv('github_ANNOVAR_hg38.hg38_multianno.txt',sep='\t')
S2=pd.read_csv('github_VEP_hg38.vcf',sep='\t', skiprows=3417)


S3=['QUAL', 'FILTER','INFO']
#'''
for f in S3:
    S1[f]='.'
#'''
for i in range(S1.shape[0]):
    for j in range(S2.shape[0]):
        if S1.Start.values[i]==S2.POS.values[j]:
           for f3 in S3:
            S1[f3].values[i]=S2[f3].values[j]
S1=S1.drop(['End'],axis=1)      

S1.to_csv('output2.csv')



