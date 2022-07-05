# PSA ~ Affine gap penalty ~ Divide conquer method

h = 2
g = 1
m = 1
mis = -2

# mathch or mismatch
def p(si, tj):
    if si == tj:
        return m
    else:
        return mis

class PSA_AGP_DC(object):
    def __init__(self, A, B):
        self.A = A
        self.B = B
        self.score = 0
        # record the match location in t or gap
        self.Align_s = [''] * len(A)
        self.traceback()

    def BestScore(self, x, s, t):
        # compute the sp score of two seqs  计算记录矩阵
        m = len(s)
        n = len(t)
        a = [0] * (n + 1)
        b = [0] * (n + 1)
        c = [0] * (n + 1)

        for j in range(1, n + 1):
            a[j] = -float('Inf')
            b[j] = -h - j * g
            c[j] = -float('Inf')
        a[0] = 0
        b[0] = -h

        for i in range(1, m + 1):
            olda = a[0]
            oldb = b[0]
            oldc = c[0]
            a[0] = -float('Inf')
            b[0] = -float('Inf')
            c[0] = -h - i * g
            for j in range(1, n + 1):
                tempa = a[j]
                tempb = b[j]
                tempc = c[j]
                if j == 1 and x != 0:
                    if x == 1:
                        if i > 1:
                            a[1] = p(s[i - 1], t[j - 1]) - h - g * (i - 1)
                        else:
                            a[1] = p(s[i - 1], t[j - 1])
                        b[1] = -float('Inf')
                        c[1] = -float('Inf')
                    elif x == 2:
                        a[1] = -float('Inf')
                        if i > 1:
                            b[1] = -h - g * i - h
                        else:
                            b[1] = -h - g
                        c[1] = -float('Inf')
                    elif x == 3:
                        a[1] = -float('Inf')
                        b[1] = -float('Inf')
                        c[1] = -h - g * i
                else:
                    c[j] = max(-h - g + a[j], -h - g + b[j], -g + c[j])
                    a[j] = p(s[i - 1], t[j - 1]) + max(olda, oldb, oldc)
                    b[j] = max(-h - g + a[j - 1], -g + b[j - 1], -h - g + c[j - 1])

                olda = tempa
                oldb = tempb
                oldc = tempc
        return a, b, c

    # Divide and conquer strategy
    def Align(self,x1, x2, a, b, c, d):
        if self.A[a:b] == '' or self.B[c:d] == '':
            if self.A[a:b]:
                self.Align_s[a:b] = ['-'] * len(self.A[a:b])
            return 0
        elif '' not in self.Align_s:
            return 0
        else:
            i = (a + b) // 2  # //取整除
            # dimension: d-c
            pref_a, pref_b, pref_c = self.BestScore(x1, self.A[a:i], self.B[c:d])
            t_rev = self.B[c:d][::-1]
            s_rev = self.A[i + 1:b][::-1]
            suff_a, suff_b, suff_c = self.BestScore(x2, s_rev, t_rev)
            suff_a = suff_a[::-1]
            suff_b = suff_b[::-1]
            suff_c = suff_c[::-1]

            posmax = c - 1  # 定位
            typemax = '-'  # 标志
            Vmax = -float('Inf')  # 值

            # find the best match (i,j) i from s, j from t   使用以上标志
            for j in range(c, d):
                maxpref = max(pref_a[j - c], pref_b[j - c], pref_c[j - c])
                maxsuff = max(suff_a[j + 1 - c], suff_b[j + 1 - c], suff_c[j + 1 - c])
                if maxpref + p(A[i], B[j]) + maxsuff > Vmax:
                    posmax = j
                    typemax = 'N'
                    Vmax = maxpref + p(A[i], B[j]) + maxsuff
                if pref_a[j - c + 1] - h - g + suff_a[j - c + 1] > Vmax:  # a - a
                    posmax = j
                    typemax = '1'
                    Vmax = pref_a[j - c + 1] - h - g + suff_a[j - c + 1]
                if pref_b[j - c + 1] - h - g + suff_b[j - c + 1] > Vmax:  # b - b
                    posmax = j
                    typemax = '2'
                    Vmax = pref_b[j - c + 1] - h - g + suff_b[j - c + 1]
                if pref_c[j - c + 1] + h - g + suff_c[j - c + 1] > Vmax:  # c - c
                    posmax = j
                    typemax = '3'
                    Vmax = pref_c[j - c + 1] + h - g + suff_c[j - c + 1]
                if pref_a[j - c + 1] - h - g + suff_b[j - c + 1] > Vmax:  # a - b
                    posmax = j
                    typemax = '4'
                    Vmax = pref_a[j - c + 1] - h - g + suff_b[j - c + 1]
                if pref_b[j - c + 1] - h - g + suff_a[j - c + 1] > Vmax:  # b - a
                    posmax = j
                    typemax = '5'
                    Vmax = pref_b[j - c + 1] - h - g + suff_a[j - c + 1]
                if pref_a[j - c + 1] - g + suff_c[j - c + 1] > Vmax:  # a - c
                    posmax = j
                    typemax = '6'
                    Vmax = pref_a[j - c + 1] - g + suff_c[j - c + 1]
                if pref_c[j - c + 1] - g + suff_a[j - c + 1] > Vmax:  # c - a
                    posmax = j
                    typemax = '7'
                    Vmax = pref_c[j - c + 1] - g + suff_a[j - c + 1]
                if pref_b[j - c + 1] - g + suff_c[j - c + 1] > Vmax:  # b - c
                    posmax = j
                    typemax = '8'
                    Vmax = pref_b[j - c + 1] - g + suff_c[j - c + 1]
                if pref_c[j - c + 1] - g + suff_b[j - c + 1] > Vmax:  # c - b
                    posmax = j
                    typemax = '9'
                    Vmax = pref_c[j - c + 1] - g + suff_b[j - c + 1]
            # i match gap
            if typemax == '1':  # a - a
                self.Align_s[i] = '-'
                self.Align(0, 1, a, i, c, posmax + 1)
                self.Align(1, 0, i + 1, b, posmax + 1, d)
            elif typemax == '2':  # b - b
                self.Align_s[i] = '-'
                self.Align(0, 2, a, i, c, posmax + 1)
                self.Align(2, 0, i + 1, b, posmax + 1, d)
            elif typemax == '3':  # c - c
                self.Align_s[i] = '-'
                self.Align(0, 3, a, i, c, posmax + 1)
                self.Align(3, 0, i + 1, b, posmax + 1, d)
            elif typemax == '4':  # a - b
                self.Align_s[i] = '-'
                self.Align(0, 1, a, i, c, posmax + 1)
                self.Align(2, 0, i + 1, b, posmax + 1, d)
            elif typemax == '5':  # b - a
                self.Align_s[i] = '-'
                self.Align(0, 2, a, i, c, posmax + 1)
                self.Align(1, 0, i + 1, b, posmax + 1, d)
            elif typemax == '6':  # a - c
                self.Align_s[i] = '-'
                self.Align(0, 1, a, i, c, posmax + 1)
                self.Align(3, 0, i + 1, b, posmax + 1, d)
            elif typemax == '7':  # c - a
                self.Align_s[i] = '-'
                self.Align(0, 3, a, i, c, posmax + 1)
                self.Align(1, 0, i + 1, b, posmax + 1, d)
            elif typemax == '8':  # b - c
                self.Align_s[i] = '-'
                self.Align(0, 2, a, i, c, posmax + 1)
                self.Align(3, 0, i + 1, b, posmax + 1, d)
            elif typemax == '9':  # c - b
                self.Align_s[i] = '-'
                self.Align(0, 3, a, i, c, posmax + 1)
                self.Align(2, 0, i + 1, b, posmax + 1, d)
            # i match j
            elif typemax == 'N':
                self.Align_s[i] = posmax
                self.Align(0, 0, a, i, c, posmax)
                self.Align(0, 0, i + 1, b, posmax + 1, d)
            return Vmax

    def traceback(self):
        #  compute the record matrix
        self.score = self.Align(0, 0, 0, len(self.A), 0, len(self.B))
        idxB_old = -1
        self.alA = ""
        self.alB = ""
        for i in range(len(self.Align_s)):
            if self.Align_s[i] == '-':
                self.alA += self.A[i]
                self.alB += '-'
            else:
                for j in range(idxB_old + 1, self.Align_s[i]):
                    self.alA += '-'
                    self.alB += self.B[j]
                self.alA += self.A[i]
                self.alB += self.B[self.Align_s[i]]
                idxB_old = self.Align_s[i]
        for i in range(idxB_old + 1, len(self.B)):
            self.alA += '-'
            self.alB += self.B[i]

def input_seq(f):
    X = []
    X_name = []
    for line in f:
        if not line.startswith('>'):
            X.append(line.replace('\n', ''))  # 去掉行尾的换行符
        else:
            X_name.append(line.replace('\n', ''))
    f.close()
    return ''.join(X), ''.join(X_name)

def output_seq(X, X_name, f):
    f.write(X_name +'\n')
    n = 0
    seq = []
    for i in X:
        seq.append(i)
        n += 1
        if n > 80:
            f.write(''.join(seq)+'\n')
            n = 0
            seq = []
    f.close()

if __name__ == "__main__":
    # input sequences
    f1 = open('C:/Users/LYZ/Desktop/1.fasta','r')
    f2 = open('C:/Users/LYZ/Desktop/2.fasta','r')
    A, A_name = input_seq(f1)
    B, B_name = input_seq(f2)

    print('序列A的长度为：', len(A))
    print('序列B的长度为：', len(B))
    psa = PSA_AGP_DC(A, B)

    # output sequences
    f3 = open('C:/Users/LYZ/Desktop/9.fasta','w')
    f4 = open('C:/Users/LYZ/Desktop/10.fasta','w')
    output_seq(psa.alA, A_name, f3)
    output_seq(psa.alB, B_name, f4)
    print('得分为：', psa.score)
    print('比对结束！\n')
