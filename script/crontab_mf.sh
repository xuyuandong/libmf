#!/bin/bash

cd $(dirname `ls -ls $0 | awk '{print $NF;}'`)/..
WK_DIR=`pwd`
cd script/

date=`date +"%Y%m%d"`

ALARM_SCRIPT=$WK_DIR/script/alarm.sh
source $WK_DIR/conf/default.conf
source $WK_DIR/script/ad_algo_fun_trap.sh

LOG_FILE=$WK_DIR/log/`date +"%Y%m%d"`.log
[ "$open_log" == "true"  ] && exec &>$LOG_FILE

set -o pipefail
set -o errexit
trap 'traperror $? $LINENO $BASH_LINENO "$BASH_COMMAND" $(printf "::%s" ${FUNCNAME[@]}) "$ALARM_SCRIPT" '  ERR   

find $WK_DIR/log/*.log -mtime +10 -delete || true
##########################################################
PYTHON="/usr/local/bin/python"

function matrix_factorization()
{
trap 'traperror $? $LINENO $BASH_LINENO "$BASH_COMMAND" $(printf "::%s" ${FUNCNAME[@]}) "$ALARM_SCRIPT" '  ERR   
app=$WK_DIR/bin/mf
input=$WK_DIR/data/$1
output=$WK_DIR/result/$2
cat $input | $app --ipopt=$WK_DIR/conf/ipopt.opt --v=2 >$output
}


timestamp=`date`
echo "begin at: $timestamp"

matrix_factorization "cust" "cust_alpha_beta"  

timestamp=`date`
echo "end at: $timestamp"

echo "successfully done $timestamp" | $PYTHON $WK_DIR/script/mailbox.py
