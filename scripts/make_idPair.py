import sys

# the input is from stdin, which is like:
# chr1	HG02587-hap1-0000005
# chr1	HG02587-hap1-0000017
# chr2	HG02587-hap1-0000028
# chr3	HG02587-hap1-0000027
# chr4	HG02587-hap1-0000004
# chr4	HG02587-hap1-0000019
# chr5	HG02587-hap1-0000007 

# the output is like:
# chr1	HG02587-hap1-0000005,HG02587-hap1-0000017
# chr2	HG02587-hap1-0000028
# chr3	HG02587-hap1-0000027
# chr4	HG02587-hap1-0000004,HG02587-hap1-0000019   

result = {}
for line in sys.stdin:
    tmp = line.strip().split()
    if len(tmp) < 2:continue
    ref, qry = tmp[0], tmp[1]
    if ref not in result:
        result[ref] = [qry,]
    else:
        result[ref].append(qry)

# 打印输出
for key, values in result.items():
    print(f"{key}\t{','.join(map(str, values))}")

