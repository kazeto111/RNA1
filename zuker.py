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

a= 6
b= -1
c= 0.1
hairpin_energy = np.array((np.inf, np.inf, 4.4, 4.6, 4.7, 4.4, 5.0, 4.5, 5.4, 5.5, 5.6,
                         5.7, 5.8, 5.9, 5.9, 6.0, 6.1, 6.1, 6.2, 6.2, 6.3, 6.3,
                         6.4, 6.4, 6.5, 6.5, 6.5, 6.6, 6.6, 6.7, np.inf), dtype=float)
internal_energy = np.array((np.inf, 1.0, 1.0, 1.1, 2.0, 2.0, 2.1, 2.3, 2.4, 2.5, 2.6, 2.7,
                 2.8, 2.9, 2.9, 3.0, 3.1, 3.1, 3.2, 3.3, 3.3, 3.4, 3.4, 3.5,
                 3.5, 3.5, 3.6, 3.6, 3.7, 3.7, np.inf), dtype=float)
bulge_energy = np.array((1.0, 1.0, 1.0, 1.1, 2.0, 2.0, 2.1, 2.3, 2.4,
                         2.5, 2.6, 2.7, 2.8, 2.9, 2.9,3.0, 3.1, 3.1, 3.2, 3.3,
                         3.3, 3.4, 3.4, 3.5, 3.5, 3.5,
                         3.6, 3.6, 3.7, 3.7, np.inf), dtype=float)

#初期化
mdp = np.full((sampleRNA_size, sampleRNA_size), np.inf, dtype=float)
vdp = np.full((sampleRNA_size, sampleRNA_size), np.inf, dtype=float)
wdp = np.zeros((sampleRNA_size, sampleRNA_size), dtype=float)

#i番目とj番目が塩基対を形成するか判断する関数
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

#stackingのenergyを計算するdp
def calc_stacking(a,b,c,d):
    if (a == 1 or b == 1) and (c == 1 or d == 1): #塩基対にUが含まれるかどうかで判断
        return -0.5
    elif (a == 1 or b == 1) or (c == 1 or d == 1):
        return -2.0
    else:
        return -3.0

#トレースバックのためのフラグ行列をつくsる
vdp_flag_matrix = np.full((sampleRNA_size, sampleRNA_size), 100, dtype=int)
for j in range(0, sampleRNA_size):
    for i in range(j, -1, -1):
        if j-i > 2 and judge_delta(sampleRNA_int[i], sampleRNA_int[j]):
            f2list = np.full((j-i-2, j-i-2), np.inf, dtype=float)
            for h in range(i+1, j-1):
                for l in range(max(h + 1, -i+j+h-32), j):
                    if judge_delta(sampleRNA_int[h], sampleRNA_int[l]):
                        if h == i+1 and l == j-1:
                            f2list[h - i - 1][l - max(h + 1, -i + j + h - 32)] \
                                = calc_stacking(sampleRNA_int[i], sampleRNA_int[j],
                                                sampleRNA_int[h], sampleRNA_int[l])

                        elif h == i+1 or l == j-1:
                            f2list[h - i - 1][l - max(h + 1, -i + j + h - 32)] \
                                = bulge_energy[h-i+j-l-3] + vdp[h][l]
                        else:
                            f2list[h-i-1][l-max(h+1, -i+j+h-32)] = \
                                internal_energy[h-i+j-l-3] + vdp[h][l]

            klist = np.full(j-i-2, np.inf, dtype=float)
            for k in range(i+1, j-1):
                klist[k-i-1] = mdp[i+1][k] + mdp[k+1][j-1]

            minlist = np.full(3, np.inf, dtype=float)
            if j-i-2 <= 30:
                minlist[0] = hairpin_energy[j-i-2]
            minlist[1] = np.amin(f2list)
            minlist[2] = np.amin(klist) + a + b
            vdp[i][j] = np.min(minlist)
            vdp_flag_matrix[i][j] = np.argmin(minlist)


            #以下mdpの再帰
            klist = np.full(j - i, np.inf, dtype=float)
            for k in range(i,j):
                klist[k-i] = mdp[i][k] + mdp[k+1][j]

            mdp[i][j] = min(vdp[i][j]+b, mdp[i+1][j]+c, mdp[i][j-1]+c, np.amin(klist))

        #以下wdpの再帰
        klist = np.full(j - i, np.inf, dtype=float)
        if i == j:
            wdp[i][j] = np.inf
        else:
            for k in range(i, j):
                klist[k - i] = wdp[i][k] + wdp[k + 1][j]
            wdp[i][j] = min(wdp[i + 1][j], wdp[i][j - 1], vdp[i][j], np.amin(klist))
            # print("wdp")
            # print("wdp[i + 1][j]", wdp[i + 1][j])
            # print("wdp[i][j - 1]", wdp[i][j - 1])
            # print("vdp[i][j]", vdp[i][j])
            # print("np.amin(klist)", np.amin(klist))

# print("wdp")
# for row in wdp:
#     print(row)
#
# print("vdp")
# for row in vdp:
#     print(vdp)
#
# print("mdp")
# for row in mdp:
#     print(row)

#以下トレースバック
#トレースバック用の関数を定義
#vdpのトレースバック用の関数
def vdp_trace(record, i, j):
    stack = []
    stack.append([i, j])
    while stack:
        i, j = stack.pop()
        if vdp_flag_matrix[i][j] == 0:
            pass
        elif vdp_flag_matrix[i][j] == 1:
            for h in range(i + 1, j - 1):
                for l in range(max(h + 1, -i + j + h - 32), j):
                    if judge_delta(sampleRNA_int[h], sampleRNA_int[l]):
                        if h == i + 1 and l == j - 1 and calc_stacking(sampleRNA_int[i], sampleRNA_int[j],
                                                                       sampleRNA_int[h], sampleRNA_int[l]) + vdp[h][
                            l] == vdp[i][j]:
                            record.append([h, l])
                            stack.append([h, l])

                        elif (h == i + 1 or l == j - 1) and bulge_energy[h - i + j - l - 3] + vdp[h][l] \
                                == vdp[i][j]:
                            record.append([h, l])
                            stack.append([h, l])

                        elif internal_energy[h - i + j - l - 3] + vdp[h][l] == vdp[i][j]:
                            record.append([h, l])
                            stack.append([h, l])
        else:
            for k in range(i + 1, j - 1):
                if mdp[i + 1][k] + mdp[k + 1][j - 1] + a + b == vdp[i][j]:
                    record.append([i+1, k])
                    record.append([k+1, j-1])
                    record = mdp_trace(i,j,k,record)
                    break
    return record


#mdpのトレースバック用の関数
def mdp_trace(i, j, k, record):
    stack = []
    stack.append([i+1, k])
    stack.append([k+1, j-1])
    while stack:
        i, j = stack.pop()
        if i>= j:
            pass
        else:
            if mdp[i][j] == vdp[i][j] + b:
                record = vdp_trace(record, i, j)
            elif mdp[i][j] == mdp[i+1][j] + c:
                stack.append([i+1, j])
            elif mdp[i][j] == mdp[i][j-1] + c:
                stack.append([i, j-1])
            else:
                for k in range(i, j):
                    # print("e")
                    if mdp[i][k] + mdp[k + 1][j] == mdp[i][j]:
                        stack.append([k + 1, j])
                        stack.append([i, k])
                        record.append([i, k])
                        record.append([k+1, j])
                        # print("f")
                        break
    return record

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
    elif wdp[i+1][j] == wdp[i][j]:
        stack.append([i+1,j])
        #print("b")
    elif wdp[i][j-1] == wdp[i][j]:
        stack.append([i,j-1])
        #print("c")
    elif judge_delta(sampleRNA_int[i],sampleRNA_int[j]) and j-i-1 >= 3 and wdp[i][j] == vdp[i][j]:
        record.append([i,j])
        newstack = []
        newstack.append([i,j])
        record = vdp_trace(record, i, j)

    else:
        for k in range(i, j-1):
            #print("e")
            if wdp[i][k] + wdp[k+1][j] == wdp[i][j]:
                stack.append([k+1,j])
                stack.append([i,k])
                #print("f")
                break
    #print(stack)


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
#print("score", wdp[0][sampleRNA_size-1])

with open("/Users/futo/Desktop/浅井先生実験/名称未設定フォルダ/kadai5_output.txt", "w") as f:
    f.write("最適構造")
    f.write(output)
    f.write("最小自由エネルギー")
    f.write(str(wdp[0][sampleRNA_size-1]))

















