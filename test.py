import re

with open('/home/dlbp/annotation-pipeline/s9_GATK_vcf/chromosome_code.txt', 'r') as ref:
    for line in ref:
        match = re.match(r'>(\w*.\d*)', line)
        with open('/home/dlbp/annotation-pipeline/s9_GATK_vcf/chromosome10_code.txt', 'a') as chro:
            if match.group(1):
                print(match.group(1) + '\n')
                chro.write(match.group(1)+'\n')
            else:
                print('not working')
