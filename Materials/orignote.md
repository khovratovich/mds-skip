## Trail of width 2: Vector

Suppose there is a 2-round differential trail:
\begin{align}
(\mathbf{A}_1,\mathbf{A}_2,0,\ldots,0) \xrightarrow{\mathcal{L}}(\mathbf{b}_1,\mathbf{b}_2,0,\ldots,0) \xrightarrow{\mathcal{S}}\\
(\mathbf{C}_1,\mathbf{C}_2,0,\ldots,0)\xrightarrow{\mathcal{L}}(\mathbf{d_1},\mathbf{d_2},0,\ldots,0)\\
\end{align}
where $\mathbf{A}$ are matrices and $\mathbf{b}$ are vectors.

Here by $\mathbf{A}\xrightarrow{\mathcal{L}}\mathbf{b}$ we mean that
$$
\mathcal{L}(\mathbf{A}\cdot\mathbf{X}) = \mathbf{b}\cdot\mathbf{X} = \sum_i \mathbf{b}[i]\cdot\mathbf{X}[i]
$$
where $\mathbf{A}\cdot\mathbf{X}$ is the matrix-vector multiplication

If the following trail condition holds:
$$
\frac{\mathbf{b_2}^r}{\mathbf{b_1}^r}\circ\mathbf{C}_1 = \mathbf{C}_2
$$
And the following condition holds for the input to the $\mathcal{S}$ layer
\begin{align}
W_1^B[0]:W_1^B[1]=\mathbf{b}_1:\mathbf{b}_2
\end{align}
Then the structure of states propagates as follows:
\begin{align}
&(W_1^A[0]+\mathbf{A}_1\cdot\mathbf{X},W_1^A[1]+\mathbf{A}_2\cdot \mathbf{X},W_1^A[2],W_1^A[3],\ldots)\xrightarrow{\mathcal{L}}\\
&(W_1^B[0]+\mathbf{b}_1\circ\mathbf{X},\underbrace{W_1^B[1]+\mathbf{b}_2\circ\mathbf{X}}_{=\frac{\mathbf{b_2}}{\mathbf{b_1}}(W_1^B[0]+\mathbf{b}_1\circ\mathbf{X})},W_1^B[2],W_1^B[3],\ldots)\xrightarrow{\mathcal{S}}\\
&((W_1^B[0]+\mathbf{b}_1\circ\mathbf{X})^r,\frac{\mathbf{b_2}^r}{\mathbf{b_1}^r}(W_1^B[0]+\mathbf{b}_1\circ\mathbf{X})^r,W_1^C[2],W_1^C[3],\ldots)\sim\\
&(W_1^C[0]+\mathbf{C}_1\cdot \mathbf{Y},\frac{\mathbf{b_2}^r}{\mathbf{b_1}^r}(W_1^C[0]+\mathbf{C}_1\cdot\mathbf{Y}),W_1^C[2],W_1^C[3],\ldots)=\\
&(W_1^C[0]+\mathbf{C}_1\cdot\mathbf{Y},W_1^C[1]+\underbrace{\frac{\mathbf{b_2}^r}{\mathbf{b_1}^r}\circ\mathbf{C}_1}_{=\mathbf{C}_2}Y,W_1^C[2],W_1^C[3],\ldots)\xrightarrow{\mathcal{L}}\\
&(W_2^B[0]+\mathbf{d_1}\circ \mathbf{Y},W_2^B[1]+\mathbf{d_2}\circ \mathbf{Y},W_2^B[2],W_2^B[3],\ldots)
\end{align}
