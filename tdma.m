% 初期化
clear all;

% Au=d
A = [3, 1, 0; 1, 4, 2; 0, 2, 5];
d = [5, 15, 19];

% 行列Aが三重対角行列か確認
Adim = size(A, 1);
n = 0;
for i = 2 : Adim - 1
    Adiags_p = spdiags(A, i);% 非三重対角部分を抽出
    Adiags_m = spdiags(A, -i);% 非三重対角部分を抽出
    n = n + nnz(Adiags_p) + nnz(Adiags_m);% 非三重対角部分の非ゼロ要素の数を数える。
end
if  n > 0
    disp('行列Aは三重対角行列ではありません。計算を終了します。')
    return;% スクリプト終了
end

% 成分抽出
b = spdiags(A, 0);% 対角成分の抽出
a = spdiags(A, - 1);% 三重対角成分の抽出（マイナス側）
c = spdiags(A, 1);% 三重対角成分の抽出（プラス側）
a = a(1 : end - 1);% 余計な成分を削除
c = c(2 : end);% 余計な成分を削除

% 求解
u = TDMA(a, b, c, d);
disp(u)

%% 以下関数

function[u] =  TDMA(a, b, c, d)

%a,cはn-1個の要素
%b,dはn個の要素

n = size(d, 2);% 未知変数の数
e = zeros(n - 1);
f = zeros(n - 1);

% 1番目の係数を求める
e(1) = c(1) / b(1);
f(1) = d(1) / b(1);

% 1からn-1番目の係数を求める
for i = 2 : n - 1
    e(i) = c(i) / (b(i) - a(i - 1) * e(i - 1));
    f(i) = (d(i) - a(i - 1) * f(i - 1)) / (b(i) - a(i - 1) * e(i - 1));
end

% n番目の解を求める
u(n) = (d(n) - a(n - 1) * f(n - 1)) / (b(n) - a(n - 1) * e(n - 1));

% n-1から1番目の解を求める
for i = n - 1 : - 1 : 1
    u(i) = f(i) - e(i) * u(i + 1);
end

end