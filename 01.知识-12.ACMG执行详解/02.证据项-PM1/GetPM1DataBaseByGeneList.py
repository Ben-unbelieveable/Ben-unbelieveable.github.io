#!/usr/bin/python
# -*- encoding: utf-8 -*-
'''
@File    :   GetPM1DataBaseByGeneList.py
@Time    :   2024/11/26 15:01:54
@Author  :   Liu.Bo
@Version :   1.0.0.0
@Contact :   liubo4@genomics.cn/614347533@qq.com
@WebSite :   http://www.ben-air.cn/
@Description :
@Notes :
@References: :
'''
from argparse import ArgumentParser
from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport
import pandas as pd
import os 
from tqdm import tqdm
import sys
import time
program = 'GetPM1DataBaseByGeneList.py'
version = '1.0.0.0'

def gnomAD_zscore_oe(gene,reference_genome="GRCh37"):
    transport = AIOHTTPTransport(url="https://gnomad.broadinstitute.org/api")
    client = Client(transport=transport, fetch_schema_from_transport=True)
    # For brevity, and to keep the focus on the Python code, we don't include every
    # field from the raw query here.
    query_seq = (
        f"""
    {{
        gene(gene_symbol: "{gene}", reference_genome: GRCh37) {{
            gene_id
            symbol
            ncbi_id
            chrom
    		start
    		stop
            gnomad_constraint {{
                oe_mis
                oe_mis_lower
                oe_mis_upper
                mis_z
            }}
            gnomad_v2_regional_missense_constraint {{
                regions {{
                    chrom
                    start
                    stop
                    obs_exp
                }}
            }}
        }}
    }}
    """
    )
    query = gql(query_seq)
    gnomad_data = client.execute(query)
    gnomad_data["input_symbol"] = gene
    return gnomad_data


def parse_info(gnomdb_info):
    df_out=pd.DataFrame(columns=["chrom", "start", "stop","ncbi_id","symbol","input_symbol","obs_exp","z-score"])
    gene_ncbi_id = gnomdb_info["gene"]["ncbi_id"]
    gene_symbol = gnomdb_info["gene"]["symbol"]
    input_symbol = gnomdb_info["input_symbol"]
    if "gene" in gnomdb_info and "gnomad_constraint" in gnomdb_info:
        z_score = gnomdb_info["gene"]["gnomad_constraint"]["mis_z"]
        regions = gnomdb_info["gene"]["gnomad_v2_regional_missense_constraint"]["regions"]
        if len(regions) > 0:
            for region in regions:
                df_out.loc[len(df_out)+1] = [region["chrom"],region["start"],region["stop"],gene_ncbi_id,gene_symbol,input_symbol,region["obs_exp"],z_score]
        else:
            df_out.loc[len(df_out)+1] = [gnomdb_info["gene"]["chrom"],gnomdb_info["gene"]["start"],gnomdb_info["gene"]["stop"],gene_ncbi_id,gene_symbol,input_symbol, gnomdb_info["gene"]["gnomad_constraint"]["oe_mis"],z_score]
        return df_out
    else:
        with open("ignore_gene.list","a") as f:
            f.write(input_symbol+"\tNo mis_z\n")
        return df_out

def deal_single_gene(gene):
    return parse_info(gnomAD_zscore_oe(gene))

def get_gene_list(file_name,colname):
    return set(pd.read_csv(file_name,sep="\t")[colname])

def save_df(df,file_name):
    if os.path.exists(file_name):
        header = False
        mode = 'a'      # 追加模式
    else:
        header = True
        mode = 'w'      # 写入模式
    df.to_csv(file_name, sep="\t", index=False, header=header, mode=mode)
def deal_list_gene(gene_list,file_name):
    total=len(gene_list)
    deal_num = 0
    success_num=0
    fail_num = 0
    #for gene in tqdm(gene_list):
    for gene in gene_list:
        deal_num+=1
        #print(str(total)+"\t:"+gene)
        try:
            parse_df = deal_single_gene(gene)
            #print(parse_df)
            save_df(parse_df,file_name)
            success_num+=1
        except:
            #print("Error: "+gene)
            with open("error_gene.txt","a") as f:
                f.write(gene+"\n")
            fail_num+=1
        info = f"deal gene num ({deal_num}/{total}); success gene num({success_num}); failed gene num({fail_num}) \tcurrent gene:\t{gene}" +" "*20
        # print(info)
        # sys.stdout.write("\r"+str(info))
        # sys.stdout.flush()
        print('\r' + info, end='', flush=True)
    return 0


def main():
    parser = ArgumentParser(prog=program)
    parser.add_argument('-in' , dest='input' , required=True, action='store', type=str, help='input database File contain gene symbol list')
    parser.add_argument('-col' , dest='colname' , default="gene_name", action='store', type=str, help='column name of gene symbol in list File')
    parser.add_argument('-i' , dest='ignore' , default="ignore_gene.list", action='store', type=str, help='ignore gene symbol list (have no oe/z-score in gnomAD)')
    parser.add_argument('-out', dest='output', required=True, action='store', type=str, help='output database File in bed format')
    args = parser.parse_args()
    GeneList_ignore = get_gene_list(args.ignore,"gene_name")
    if os.path.exists(args.output):
        GeneList_inInput = get_gene_list(args.input,args.colname)
        GeneList_inOutput = get_gene_list(args.output,"input_symbol")
        GeneList_toDeal = GeneList_inInput-GeneList_inOutput-GeneList_ignore
    else:
        GeneList_toDeal = get_gene_list(args.input,args.colname)
    print("Total gene num(need to run) is "+str(len(GeneList_toDeal)) +  " ; Done genenum:"+str(len(GeneList_inOutput))+ " ; ignore genenum:"+str(len(GeneList_ignore)))
    deal_list_gene(GeneList_toDeal,args.output)
    Miss_num_inlastCircle = len(GeneList_toDeal)
    circle=1
    while Miss_num_inlastCircle>0:
        time.sleep(10)
        circle = circle+1
        Miss_gene = GeneList_inInput-get_gene_list(args.output,"input_symbol")-GeneList_ignore
        Miss_num=len(Miss_gene)
        if len(Miss_gene)>0:
            deal_list_gene(GeneList_toDeal,args.output)
        Miss_num_inlastCircle = Miss_num
        print("\nMiss gene num after circle"+str(circle)+" is "+str(Miss_num_inlastCircle))


if __name__ == '__main__':
    main()

