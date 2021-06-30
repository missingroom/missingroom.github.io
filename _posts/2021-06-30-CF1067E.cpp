---
title: CF1067E Random Forest Rank
tags: 线性代数 矩阵 树 dp 概率期望
---

考虑矩阵的秩的意义，即非零子式的最大阶数，也就是说，我们要找到一个最大的子式使得其行列式不为零。由于矩阵是无向图的邻接矩阵，若子式行列式不为零，选择包括选择的所有行列下标的主子式也一定行列式不为零，而这个主子式阶数更大，所以我们只用考虑主子式。邻接矩阵的主子式既是其导出子图的邻接矩阵，所以问题变为什么样的森林行列式不为零。

考虑行列式的定义：
$$
\det A=\sum\limits_{\pi}(-1)^{\tau(\pi)}\prod\limits_{i=1}^nA_{i,\pi_i}
$$
注意到如果 $A_{i,\pi_i}=1$，则存在边 $(i,\pi_i)$，由于无自环，该排列不存在长度为 $1$ 的循环，由于无环，该排列不存在长度大于 $2$ 的循环，所以排列只存在长度为 $2$ 的循环，也就意味着，当且仅当 $\pi$ 为这张图的完美匹配时，$\prod\limits_{i=1}^nA_{i,\pi_i}\ne 0$。 

森林有一个显然的性质，完美匹配不多于一种，考虑从叶子开始匹配，匹配方法显然唯一。

所以森林的邻接矩阵行列式不为零，当且仅当存在完美匹配；一个森林的邻接矩阵的秩即其最大匹配的二倍。

问题转变为求这个森林最大匹配的期望。

设 $f_i$ 为 $i$ 与其儿子匹配的方案数，答案即为 $2\sum\limits_{i=1}^nf_i$，具体的转移可以看代码，代码里的 `f[i]` 是 $\frac{f_i}{2^{n-sz_i}}$。

```cpp
#include<cstdio>
#include<vector>
int const mod=998244353,inv2=(mod+1)/2;
std::vector<int> v[500010];
int n,f[500010],pw[500010],ans,sz[500010];
void dfs(int x,int fa){
	int p=1;
	sz[x]=1;
	for(auto u:v[x])if(u!=fa){
		dfs(u,x);
		p=1ll*p*(pw[sz[u]-1]+f[u])%mod;
		sz[x]+=sz[u];
	}
	ans=(ans+1ll*pw[n-sz[x]]*(f[x]=(pw[sz[x]-1]-p+mod)%mod))%mod;
}
int main(){
	scanf("%d",&n);
	pw[0]=1;
	for(int i=1;i<=n;i++) pw[i]=pw[i-1]*2%mod;
	for(int i=1,x,y;i<n;i++)scanf("%d%d",&x,&y),v[x].push_back(y),v[y].push_back(x);
	dfs(1,0);
	printf("%d\n",2*ans%mod);
	return 0;
}
```

