#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# @Time    : 2018/1/16 下午8:13
# @Author  : chenyuelong
# @Mail    : yuelong_chen@yahoo.com
# @File    : getSubClonePipe.py
# @Software: PyCharm

import os, sys
import configparser
import argparse

path = os.path.abspath(os.path.dirname(__file__))


def getConfig(configfile):
    config = configparser.ConfigParser()
    config.read(configfile)
    return config


def getArgs():
    parser = argparse.ArgumentParser('亚克隆预测pipeline')
    parser.add_argument('--config', '-config', dest='sconfig', action='store',
                        help='样本配置文件，详情见例子', required=True)
    args = parser.parse_args()
    return args


def getStates(state):
    return state.split(',')


def preOneLib(config, substate):
    snp = config.get(substate, 'snp')
    indel = config.get(substate, 'indel')
    driver = config.get(substate, 'driver')
    return snp, indel, driver


def get_sd2tsv(perl, sd2tsv, snp, indel, out):
    cmd = '{} {} {} {} {}'.format(perl, sd2tsv, snp, indel, out)
    return cmd


def get_anno2tsv(python, anno2tsv, files, samples, ponfile, out):
    cmd = '{0} {1} -infiles  {2} -s {3} --filter {4} -o {5}'.format(python, anno2tsv, files, samples, ponfile, out)
    return cmd


def get_driver(perl, scr, table, driver, out):
    cmd = '{} {} {} {} {}'.format(perl, scr, table, driver, out)
    return cmd


def get_clonevol(rscript, clovR, cloneInput, outdir, sump, alpha):
    # outdir = '{}/clonevol_{}_{}'.format(outdir, sump, alpha)
    # os.makedirs(outdir, exist_ok=True)
    cmd = '{} {} {} {} {} {}'.format(rscript, clovR, cloneInput, outdir, sump, alpha)
    return cmd


def get_pyclone(pyclone, state, outdir):
    infile = ''
    for i in state:
        infile = '{0} {1}/{2}.tsv'.format(infile, outdir, i)

    cmd = '{0} setup_analysis ' \
          '--in_files {1} ' \
          '--working_dir {2} ' \
          '--samples {3} ' \
          '--prior total_copy_number &&' \
          '{0} run_analysis ' \
          '--config_file {2}/config.yaml --seed 1&& ' \
          '{0} build_table ' \
          '--config_file {2}/config.yaml ' \
          '--out_file {2}/table.old_style ' \
          '--table_type old_style  --max_clusters 6'.format(pyclone, infile, outdir, ' '.join(state))
    return cmd


def get_pyclone_plot(pyclone, outdir):
    pydict = {}
    pydict['plot_clusters'] = ['density', 'parallel_coordinates', 'scatter']
    pydict['plot_loci'] = ['density', 'parallel_coordinates', 'scatter',
                           'similarity_matrix', 'vaf_parallel_coordinates', 'vaf_scatter']
    cmd = 'echo "start plot"'
    for i in pydict.keys():
        for j in pydict[i]:
            cmd = '{0} && {1} {2} --config_file {3}/config.yaml --plot_file {3}/{2}_{4}.pdf --plot_type {4}' \
                .format(cmd, pyclone, i, outdir, j)

    return cmd


def main():
    args = getArgs()
    sampleConfig = getConfig(args.sconfig)
    # sampleConfig = getConfig('./example_sample_config.ini')
    progrmConfig = getConfig('{}/configure.ini'.format(path))

    ##get program
    perl = progrmConfig.get('program', 'perl')
    python = progrmConfig.get('program', 'python')
    rscript = progrmConfig.get('program', 'rscript')
    pyclone = progrmConfig.get('program', 'pyclone')

    ##get script
    sd2tsv_s = progrmConfig.get('scripts', 'sd2tsv')
    clonevolR_s = progrmConfig.get('scripts', 'clonevolR')
    getDriver_s = progrmConfig.get('scripts', 'getDriver')
    anno2tsv_s = progrmConfig.get('scripts', 'anno2tsv')

    ## get database
    ponfile = progrmConfig.get('database', 'ponfile')

    ## get sample config
    samplename = sampleConfig.get('sample', 'name')
    state = getStates(sampleConfig.get('sample', 'state'))
    outdir = sampleConfig.get('workdir', 'outdir')
    os.makedirs(outdir, exist_ok=True)

    ## get cmds
    cmds = []
    driver = ''
    files = ''
    samples = ''
    for i in state:
        snp, indel, driver = preOneLib(sampleConfig, i)
        files = '{} {},{}'.format(files, snp, indel)
        samples = '{} {}'.format(samples, i)
        # cmds.append(get_sd2tsv(perl, sd2tsv_s, snp, indel, '{}/{}.tsv'.format(outdir, i)))
    cmds.append(get_anno2tsv(python, anno2tsv_s, files, samples, ponfile, outdir))
    cmds.append(get_pyclone(pyclone, state, outdir))
    cmds.append(get_driver(perl, getDriver_s, '{}/table.old_style'.format(outdir),
                           driver, '{}/cloneEvaInput.txt'.format(outdir)))
    cmds.append(get_pyclone_plot(pyclone, outdir))
    for i in range(1, 100, 5):
        for j in range(1,100,5):
            sump=i/100
            alpha =j/100
            cmds.append(get_clonevol(rscript, clonevolR_s, '{}/cloneEvaInput.txt'.format(outdir),
                                     outdir,sump,alpha))

    ## writer cmds
    for i in cmds:
        print(i)


if __name__ == '__main__':
    main()
