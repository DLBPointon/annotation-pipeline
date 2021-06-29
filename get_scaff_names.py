import re
import sys
import os

listy = []

with open('./chromosome_code.txt', 'r') as ref:
    for line in ref:
        match = re.match(r'>(\w*.\d*)', line)
        m = match.group(1)
        if m:
            listy.append(m)

with open('./chro_output.txt', 'a') as chro:
    for i in listy:
        if i != listy[-1]:
            chro.write(i + ',')
        else:
            chro.write(i)
