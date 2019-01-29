#!/bin/bash

echo "Running Integration Test: unfolded+tagged ..."

SOFOS=$1
VCF=test-default.vcf

output=$(${SOFOS} -q -a 0.1 -b 1 -n 10 -u -t ${VCF}) 

hnum=$(echo "${output}" | grep '^#' | wc -l)
if [[ 10 -ne "${hnum}" ]]; then
	echo "  ERROR: Heading has ${hnum} lines instead of 10."
	exit 1
fi
header=$(echo "${output}" | awk 'NR==11')
if [[ $header != "Number,Prior,Observed,Posterior" ]]; then
	echo "  ERROR: Header has wrong column names."
	exit 1
fi

# Check posterior
expected="0.341830374360738
0.250675607864541
0.169206035308565
0.107597683991088
0.0643344485530045
0.0357933477403989
0.018194951768036
0.00820216873352734
0.00311426094101116
0.000899675382958782
0.000151445356131395"

echo "${output}" | grep -v '^#' | awk -F, 'NR>1 {print $4}' \
	| paste -d, - <(echo "${expected}") \
	| awk -F, -v m=1e-7 'function abs(v) {return v < 0 ? -v : v} {a=abs(($1-$2))/$2; if(a > m){m = a} } END{exit (m > 1e-7) }' 
if [[ $? -ne 0 ]]; then
	echo "  ERROR: Maximum relative error in posterior exceeds 1e-7."
	exit 1
fi

# Check prior
expected="0.751621782268415
0.0751621782268415
0.0413391980247626
0.0289374386173338
0.0224265149284337
0.0183897422413156
0.0156312809051183
0.0136215447887459
0.012089121000012
0.0108802089000108
0.00990099009900991"

echo "${output}" | grep -v '^#' | awk -F, 'NR>1 {print $2}' \
	| paste -d, - <(echo "${expected}") \
	| awk -F, -v m=1e-7 'function abs(v) {return v < 0 ? -v : v} {a=abs(($1-$2))/$2; if(a > m){m = a} } END{exit (m > 1e-7) }' 
if [[ $? -ne 0 ]]; then
	echo "  ERROR: Maximum relative error in prior exceeds 1e-7."
	exit 1
fi

# Check obs
expected="0
0
1
0
0
0
0
0
0
0"

echo "${output}" | grep -v '^#' | awk -F, 'NR>1 {print $3}' \
	| paste -d, - <(echo "${expected}") \
	| awk -F, -v m=1e-7 'function abs(v) {return v < 0 ? -v : v} {a=abs(($1-$2)); if(a > m){m = a} } END{exit (m > 1e-7) }' 
if [[ $? -ne 0 ]]; then
	echo "  ERROR: Maximum error in obs exceeds 1e-7."
	exit 1
fi

# Check indexes
expected="0
1
2
3
4
5
6
7
8
9
10"

echo "${output}" | grep -v '^#' | awk -F, 'NR>1 {print $1}' \
	| paste -d, - <(echo "${expected}") \
	| awk -F, -v m=1e-7 'function abs(v) {return v < 0 ? -v : v} {a=abs(($1-$2)); if(a > m){m = a} } END{exit (m > 1e-7) }' 
if [[ $? -ne 0 ]]; then
	echo "  ERROR: Maximum error in obs exceeds 1e-7."
	exit 1
fi
