import numpy as np

fasta_symbol = input()
#sampleRNAの入力
sampleRNA = ""
while True:
    i = input()
    if len(i) == 0:
        break
    else:
        sampleRNA += i
    #print(sampleRNA)
sampleRNA_size = len(sampleRNA)
#print("sampleRNA_size",sampleRNA_size)
alphabet = ["A", "U", "G", "C"]
#alphabetのインデックスでsampleRNAを数字に変換したsample_RNA_intを作成
sampleRNA_int = np.zeros((sampleRNA_size), dtype=int)

for i in range(sampleRNA_size):
    sampleRNA_int[i] = alphabet.index(sampleRNA[i])
# print(sampleRNA_int)

#dpの初期化
nussinov_dp = np.zeros((sampleRNA_size, sampleRNA_size), dtype=int)
#deltaの計算のための関数
def judge_delta(a, b):
    if a == 0 and b == 1:
        return True
    elif a == 1 and b == 0:
        return True
    elif a == 2 and b == 3:
        return True
    elif a == 3 and b == 2:
        return True
    elif a == 1 and b == 2:
        return True
    elif a == 2 and b == 1:
        return True
    else:
        return False
#再帰
for j in range(1, sampleRNA_size):
    for i in range(j-1, -1, -1):
        klist = np.zeros(j-i, dtype=int)
        for k in range(i, j):
            klist[k-i] = nussinov_dp[i][k] + nussinov_dp[k+1][j]
        if j - i - 1 < 3:
            delta = 0
        else:
            if judge_delta(sampleRNA_int[i], sampleRNA_int[j]):
                delta = 1
            else:
                delta = 0

        nussinov_dp[i][j] = max(nussinov_dp[i+1][j], nussinov_dp[i][j-1], nussinov_dp[i+1][j-1]+delta, np.max(klist))

#print(nussinov_dp)

#トレースバック
stack = []
record = []
stack.append([0,sampleRNA_size-1])
while stack:
    i,j = stack.pop(-1)
    #print(i,j)
    #print("nussinovij", nussinov_dp[i][j])
    #print("nussinovi+1j-1", nussinov_dp[i+1][j-1])
    #print(judge_delta(sampleRNA_int[i],sampleRNA_int[j]))
    if i >= j:
        #print("a")
        pass
    elif nussinov_dp[i+1][j] == nussinov_dp[i][j]:
        stack.append([i+1,j])
        #print("b")
    elif nussinov_dp[i][j-1] == nussinov_dp[i][j]:
        stack.append([i,j-1])
        #print("c")
    elif judge_delta(sampleRNA_int[i],sampleRNA_int[j]) and j-i-1 >= 3 and nussinov_dp[i+1][j-1] + 1 == nussinov_dp[i][j]:
        record.append([i,j])
        stack.append([i+1, j-1])
        #print("d")
    else:
        for k in range(i, j-1):
            #print("e")
            if nussinov_dp[i][k] + nussinov_dp[k+1][j] == nussinov_dp[i][j]:
                stack.append([k+1,j])
                stack.append([i,k])
                #print("f")
                break
    #print(stack)

#print("record", record)

#以下構造の表記
#"("で表記する所を-1で")"で表記する所を1で表現
structure_int = np.zeros(sampleRNA_size, dtype=int)

for i in range(len(record)):
    structure_int[record[i][0]] = -1
    structure_int[record[i][1]] = 1

structure = ["a"] * sampleRNA_size
for i in range(sampleRNA_size):
    if structure_int[i] == -1:
        structure[i] = "("
    elif structure_int[i] == 0:
        structure[i] = "."
    else:
        structure[i] = ")"

output = "".join(structure)
#print(output)
with open("", "w") as f:
    f.write(output)

