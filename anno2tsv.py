#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2018/1/18 下午3:08
# @Author  : chenyuelong
# @Mail    : yuelong_chen@yahoo.com
# @File    : anno2tsv.py
# @Software: PyCharm

import os, sys
import argparse
import re


def getArgs():
    parser = argparse.ArgumentParser(prog='anno2tsv.py', description='转化multianno.txt结果到tsv')
    parser.add_argument('--input_files', '-infiles', dest='infiles',
                        nargs='+', required=True, action='store',
                        help='输入文件，一个样本snp和indel结果以"，"分割，不同样本的以" "（空格）分割')
    parser.add_argument('--samples', '-s', dest='samples',
                        nargs='+', required=True, action='store',
                        help='输入文件标记，空格分割')
    parser.add_argument('--outdir','-o',dest='outdir',required=True,
                        help='输出目录',action='store')
    parser.add_argument('--driver', '-driver', dest='driver',
                        action='store', default='',
                        help='driver gene file')
    parser.add_argument('--filter', '-f', dest='filter',
                        action='store', default='',
                        help='需要过滤位点，原始vcf格式')
    # parser.print_help()
    args = parser.parse_args()
    return args


def getLineInfo(line):
    '''
    为了方便以后修改变异筛选条件，将筛选写了一个方法
    :param line: multianno.txt中的一行
    :return: 提取的有用的信息 (绝对位置，info，ref，alt)
    '''
    cells = line.strip('\n').split('\t')
    flag = '{}_{}_{}_{}'.format(cells[55], cells[56], cells[58], cells[59])
    pattern = re.compile('[\d]/[\d]:([\d]+),([\d]+)')
    gps = pattern.findall(cells[-1])
    if len(gps) > 0:
        (ref, alt) = gps[0]
    else:
        (ref, alt) = (0, 0)
        print('-------------------' + line)
    info = '{}_{}'.format(cells[7], cells[9])
    # print(flag,info,ref,alt)
    return flag, info, ref, alt


def getPon(ponfile):
    '''
    pon文件建立dict,key为绝对位置
    :param ponfile: 自建PON文件
    :return: pondict
    '''
    pondict = {}
    if ponfile == '':
        return pondict
    with open(ponfile, 'r') as f:
        lines = f.readlines()
        for line in lines:
            cells = line.strip('\n').split('\t')
            flag = '{}_{}_{}_{}'.format(cells[0], cells[1], cells[2], cells[3])
            pondict[flag] = cells[4]
    return pondict


def getmutations(files, tags):
    '''
    整理multianno.txt结果
    :param files: args.infiles
    :param tags: args.samples
    :return: 存储结果的dict,location
    '''
    mutations = {}
    location = {}
    usetag = []
    for file, tag in zip(files, tags):
        usetag.append(tag)
        splitfiles = file.split(',')
        for sf in splitfiles:
            with open(sf, 'r') as f:
                lines = f.readlines()
                lines.pop(0)  # 取出表头
                for line in lines:
                    flag, info, ref, alt = getLineInfo(line)
                    mutationFlag = '{}_{}'.format(tag, flag)
                    mutations[mutationFlag] = '{},{}'.format(ref, alt)
                    if flag not in location.keys():
                        location[flag] = info
                    for subut in usetag:
                        tflag = '{}_{}'.format(subut, flag)
                        if tflag not in mutations.keys():
                            mutations[tflag] = '1000,1'
                    for lo in location.keys():
                        tflag = '{}_{}'.format(tag, lo)
                        if tflag not in mutations.keys():
                            mutations[tflag] = '1000,1'
    return mutations, location


''' driver file中的格式
NM_005157(<i>ABL1</i>):c.758A>T(p.Y253F)	Exon 4	43.46%	血液肿瘤综合性	xyb1ABL011
NM_005157(<i>ABL1</i>):c.944C>T(p.T315I)	Exon 6	42.24%	血液肿瘤综合性	xyb1ABL003
NM_001127208(<i>TET2</i>):c.3595-2A>G	/	62.35%	血液肿瘤综合性（111个基因）	xyb1TET001
NM_001127208(<i>TET2</i>):c.3954+1G>A	/	3.05%	血液肿瘤综合性（111个基因）	xyb1TET001
'''


def getDriver(driverfile):
    '''
    通过driverFile获取driverDict
    :param driverfile: driverFile
    :return: driverdict
    '''
    driverdict = {}
    if driverfile == '':
        return driverdict
    else:
        with open(driverfile, 'r') as f:
            lines = f.readlines()
            pat1 = re.compile('<i>(\w+)</i>.+\((p\..+)\)')
            pat2 = re.compile('<i>(\w+)</i>.+(c.[\d]+[+-][1-3].+)')
            for line in lines:
                cells = line.strip('\n').split('\t')
                gps1 = pat1.findall(cells[0])
                gps2 = pat2.findall(cells[0])
                # print('line:'+line+'|||'+cells[0])
                # print(gps1)
                # print(gps2)
                if len(gps1) > 0:
                    driverdict['{}_{}'.format(gps1[0][0], gps1[0][1])] = 1
                    # print('{}_{}'.format(gps1[0][0],gps1[0][1]))
                if len(gps2) > 0:
                    # print('{}_{}'.format(gps2[0][0],gps2[0][1]))
                    driverdict['{}_{}'.format(gps2[0][0], gps2[0][1])] = 1
                    # print('{}_splicing'.format(gps2[0][0]))
    return driverdict

def filterlocation(pondict,location):
    '''
    利用PON过滤位点
    :param pondict: pondict
    :param location:
    :return: filterlovation
    '''
    for l in location.keys():
        if l in pondict.keys():
            del location['l']
    return location

def writerfile(outdir,tag,location,mutation):
    '''
    输出文件
    :param outdir:
    :param tag:
    :param location:
    :param mutation:
    :return: nothing
    '''
    with open('{}/{}.tsv'.format(outdir,tag),'w') as f:
        title = 'mutation_id\tref_counts\tvar_counts\tnormal_cn\tminor_cn' \
                '\tmajor_cn\tvariant_case\tvariant_freq\tgenotype\n'
        f.write(title)
        for l in location:
            flag = '{}_{}'.format(tag,l)
            ref,alt = mutation[flag].split(',')
            outstr = '{0}\t{1}\t{2}\t2\t0\t2\t{3}\t{4}\t{3}\n'.\
                format(l,ref,alt,location[l],int(alt)/(int(alt)+int(ref)))
            f.write(outstr)
    print('{} finish Write!'.format(tag))

        


def mytest():
    # pat1 = re.compile('<i>(\w+)</i>.+\((p\..+)\)')
    # pat2 = re.compile('<i>(\w+)</i>.+([+-])')
    #
    # line ='NM_005157(<i>ABL1</i>):c.758A>T(p.Y253F)	Exon 4	43.46%	血液肿瘤综合性	xyb1ABL011'
    #
    # line2 = 'NM_001127208(<i>TET2</i>):c.3954+1G>A	/	3.05%	血液肿瘤综合性（111个基因）	xyb1TET001'
    # gps = pat2.findall(line)
    # print(gps)
    # gps = pat1.findall(line)
    # print(gps)
    # print(getDriver('/home/chenyl/TMPS/driver.txt'))
    args = getArgs()
    print(args.infiles)
    print(args.samples)
    tags = args.samples
    test, location = getmutations(args.infiles, args.samples)
    for l in location.keys():
        t = '{}\t{}'.format(l, location[l])
        for ta in tags:
            # print(ta)
            t = '{}\t{}'.format(t, test['{}_{}'.format(ta, l)])
        # print(t)
    for tag in tags:
        writerfile('/home/chenyl/TMPS',tag,location,test)

def run():
    args = getArgs()
    mutations,locations = getmutations(args.infiles, args.samples)
    pondict = getPon(args.filter)
    locations = filterlocation(pondict,locations)
    for tag in args.samples:
        writerfile(args.outdir,tag,locations,mutations)
    print('All finish')


def main():
    # getArgs()
    run()


if __name__ == '__main__':
    main()
