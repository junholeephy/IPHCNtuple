awk '{print $1}' 3l_SR_80X_Summer16.txt > 3l_c1LLR.txt
awk '{print $1}' sync_3l_SR.txt         > 3l_c1IPHC.txt
diff 3l_c1LLR.txt 3l_c1IPHC.txt > 3l_diff.txt

