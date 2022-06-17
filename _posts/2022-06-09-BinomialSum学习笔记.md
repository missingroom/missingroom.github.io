---
title: 载谭 Binomial Sum 学习笔记
tags: 数学
---

sto EntropyIncreaser orz

我们称一个形式幂级数 $A(x)$ 是微分有限的，当且仅当存在多项式数列 $Q_0(x),Q_1(x),\cdots,Q_m(x)$，满足 $Q_m(x)\ne 0$，且 $\sum_{i=0}^mQ_i(x)A^{(i)}(x)=0$。

## 问题描述

对于一个微分有限的生成函数 $F(x)$，和一个生成函数 $G(x)$，和一个数列 $a$，如果我们对于每个 $0\le k\le n$ 知道 $\sum_{j=0}^na_j[x^j]G(x)^k$，那么我们可以在 $\Theta(n)$ 内计算出 $\sum_{j=0}^na_j[x^j]F(G(x))$。

## 解决方法

设 $F(x)$ 对应的微分方程为 $\sum_{i=0}^kQ_i(x)F^{(i)}(x)=0$，这里我们认为 $k$ 和 $Q_i(x)$ 的次数是常数，设 $G(x)$ 的常数项为 $c$。

考虑将 $F(G(x))$ 展开成 $F(G(x)-c+c)$ 的形式，并根据 $G(x)-c$ 没有常数项，我们可以对 $F$ 进行截断，即考虑求出 $F(x+c)\bmod x^{n+1}$。 

设 $F_0(x+c)=F(x+c)\bmod x^{n+1}$。

则 $\sum_{j=0}^na_j[x^j]F_0(G(x))$ 就是答案。 

注意到，由于 $\sum_{i=0}^kQ_i(x)F^{(i)}(x)=0$，则 $\sum_{i=0}^kQ_i(x+c)F^{(i)}(x+c)=0$。

将 $F$ 换成 $F_0$ 的话，第 $n+1,n+2,\cdots,n+k$ 都变成了 $0$，等式右边需要补一个 $k$ 项的多项式来修正，即 $\sum_{i=0}^kQ_i(x+c)F_0^{(i)}(x+c)=D(x)$，也即 $\sum_{i=0}^kQ_i(x)F_0^{(i)}(x)=D(x-c)$。

这样就可以 $\Theta(n)$ 递推出 $F_0(x)$ 的每一项系数了，不妨设为 $f_0$ 到 $f_n$。

那么答案就是 $\sum_{j=0}^na_j[x^j]\sum_{i=0}^nf_iG(x)^i=\sum_{i=0}^nf_i\sum_{j=0}^na_j[x^j]G(x)^i$，我们就可以 $\Theta(n)$ 求出答案了。

注意到，如果把问题描述中的 $[x^j]$ 变成 $[\frac{x^j}{j!}]$，整个过程也是没有问题的。

## 例题

### 线性插值

对于一个 $n$ 次多项式 $f$，给出 $n+1$ 个数 $f(0),f(1),\cdots,f(n)$，和一个数 $k$，求 $f(k)$。

我们对于 $0\le i\le n$，知道 $\sum_{j=0}^nf_ji^j=\sum_{j=0}^nf_j[\frac{x^j}{j!}](e^x)^i$，想要求解 $\sum_{j=0}^nf_jk^j=\sum_{j=0}^nf_j[\frac{x^j}{j!}](e^x)^k$。

问题描述中的 $a$ 数列就是 $f$ 的系数，$G(x)$ 取 $e^x$，$F(x)$ 就是 $x^k$。

$F(x)$ 满足的微分方程是 $xF^{\prime}(x)-kF(x)=0$。

设 $F_0(x+1)=F(x+1)\bmod x^{n+1}$。
$$
(x+1)F^{\prime}(x+1)-kF(x+1)=0\\
(x+1)F_0^{\prime}(x+1)-kF_0(x+1)=(n-k)\binom knx^n\\
xF_0'(x)-kF_0(x)=(n-k)\binom kn(x-1)^n\\
$$

提取系数即可，也就是说，我们设 $F_0$ 的 $i$ 次项系数为 $f_i$，有

$$
if_i-kf_i=(n-k)\binom kn(-1)^{n-i}\binom ni
$$

### P5907

给定 $n,q,k$，求 $\sum_{i=0}^ni^kq^i$。目标复杂度 $O(k)$。

$G(x)$ 取 $qe^x$，$F(x)$ 取 $\frac{1-x^{n+1}}{1-x}$，我们已知 $\forall 0\le i\le k,[\frac{x^k}{k!}](qe^x)^i=q^ii^k$，求 $\sum_{i=0}^ni^kq^i=[\frac{x^k}{k!}]F(G(x))$。

$F(x)$ 满足的微分方程是 $F(x)+(x-1)F^{\prime}(x)=(n+1)x^{n}$。

即 $-nF(x)+((2-n)x+n)F^{\prime}(x)+(x^2-x)F^{\prime\prime}(x)=0$​。

设 $F_0(x+q)=F(x+q)\bmod x^{k+1}$。

$$
-nF(x+q)+((2-n)(x+q)+n)F^{\prime}(x+q)+((x+q)^2-(x+q))F^{\prime\prime}(x+q)=0\\
-nF(x+q)+((2-n)x+2q-nq+n)F^{\prime}(x+q)+(x^2+(2q-1)x+q^2-q)F^{\prime\prime}(x+q)=0\\
-nF_0(x+q)+((2-n)x+2q-nq+n)F_0^{\prime}(x+q)+(x^2+(2q-1)x+q^2-q)F_0^{\prime\prime}(x+q)=ax^k+bx^{k-1}\\
$$

设 $c=[x^k]F(x+q),d=[x^{k-1}]F(x+q)$，则上式 $a,b$ 分别为

$$
a=-nc+(2-n)kc+k(k-1)c\\
b=-nd+(2-n)(k-1)d+(2q-nq+n)kc+(k-1)(k-2)d+(2q-1)k(k-1)c
$$

考虑求出 $c,d$，$c=\sum_{i=0}^n[x^k] (x+q)^i=\sum_{i=k}^nq^{i-k}\binom ik$。

$$
h_k=\sum_{i=k}^nq^{i-k}\binom ik\\
=\sum_{i=k}^nq^{i-k}(\binom {i-1}{k-1}+\binom{i-1}{k})\\
=\sum_{i=k-1}^{n-1}q^{i-(k-1)}\binom {i}{k-1}+q\sum_{i=k}^{n-1}q^{i-k}\binom{i}{k}\\
=h_{k-1}-q^{n-k+1}\binom n{k-1}+q(h_k-q^{n-k}\binom nk)
$$

递推即可，注意特判 $q=1$。

然后根据

$$
-nF_0(x)+((2-n)x+n)F_0^{\prime}(x)+(x^2-x)F_0^{\prime\prime}(x)=a(x-q)^k+b(x-q)^{k-1}
$$

求出 $F_0$ 即可，设 $[x^i]F_0(x)=f_i$，那么有

$$
-nf_i+(2-n)if_{i}+n(i+1)f_{i+1}+i(i-1)f_i-(i+1)if_{i+1}=a(-q)^{k-i}\binom ki+b(-q)^{k-1-i}\binom {k-1}i\\
f_i=\frac{a(-q)^{k-i}\binom ki+b(-q)^{k-1-i}\binom {k-1}i+(i-n)(i+1)f_{i+1}}{(i-n)(i+1)}
$$

我们已知 $f_{k+1}=0$，就可以求出 $f_{0\sim k}$ 了，那么答案就是 $\sum_{i=0}^ni^kq^i=\sum_{i=0}^kf_i[\frac{x^k}{k!}](qe^x)^i=\sum_{i=0}^kf_iq^ii^k$。

注意， $n\le k$ 的时候求解 $f_i$ 时的 $(i-n)$ 会出问题，暴力即可。

### P6667

给出一个 $m$ 次多项式 $h(x)$ 在 $0\sim m$ 处的点值，和两个整数 $n,q$，求 $\sum_{k=0}^n h(k)\binom nkq^k(1-q)^{n-k}$。目标复杂度 $O(m)$。

我们已知的是 $\forall 0\le k\le m,\sum_{i=0}^m h_i k^i=\sum_{i=0}^m h_i[\frac{x^i}{i!}] (e^x)^k$。

首先我们要求的就是 $\sum_{i=0}^mh_i\sum_{k=0}^nk^i\binom nkq^k(1-q)^{n-k}=\sum_{i=0}^mh_i[\frac{x^i}{i!}]\sum_{k=0}^n\binom nk(qe^x)^k(1-q)^{n-k}$。

设 $G(x)=e^x$，$F(x)=\sum_{k=0}^n\binom nk(qx)^k(1-q)^{n-k}=(qx+1-q)^n$。

$F(x)$ 满足的微分方程是 $(qx+1-q)F^{\prime}(x)-nqF(x)=0$。

设 $F_0(x+1)=F(x+1)\bmod x^{m+1}$。

$$
(q(x+1)+1-q)F^{\prime}(x+1)-nqF(x+1)=0\\
(qx+1)F^{\prime}(x+1)-nqF(x+1)=0\\
(qx+1)F_0^{\prime}(x+1)-nqF_0(x+1)=ax^m\\
$$
设 $b=[x^m]F(x+1)=[x^m](qx+1)^n=\binom nmq^m$，那么 $a=qmb-nqb=(m-n)\binom nmq^{m+1}$。

所以我们有
$$
(qx+1-q)F_0^{\prime}(x)-nqF_0(x)=a(x-1)^m\\
$$
设 $[x^i]F_0(x)=f_i$，那么有
$$
qif_i+(1-q)(i+1)f_{i+1}-nqf_i=a(-1)^{m-i}\binom mi\\
f_i=\frac{a(-1)^{m-i}\binom mi+(q-1)(i+1)f_{i+1}}{q(i-n)}
$$
递推即可，注意， $n\le m$ 的时候求解 $f_i$ 时的 $(i-n)$ 会出问题，暴力即可。
