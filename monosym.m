%lexicographic functions, 6.45

format long

global nCk
n = 67;
gen_comb_table(n);
fprintf("table made\n")
%check_decr(2, n)
%carr = gen_monosym(10,n);
%inf = influence(carr, n);
%ent = entropy(carr,n);
%testyCoefs = make_coef_table(n)

testyC = make_C_table(n)
plot(testyC)

%check decreasing
function [] = check_decr(start, upto)
    for m = start:upto
        tab = make_C_table(m);
        val = tab(1);
        for i = 2:length(tab)
            if tab(i) >= val
                fprintf("ERROR AT level " + m + "\n")
                fprtinf(tab)
                break
            else
                val = tab(i);
            end
        end
        fprintf("good at level" + m + "\n")
    end
end

function C = calc_C(coef_arr,n)
    C = double(entropy(coef_arr,n)/influence(coef_arr,n));
end

function C_table = make_C_table(n)
    global nCk
    C_table = zeros(1, floor((n+1)/2));
    coef_table = make_coef_table(n);
    for w = 1:floor((n+1)/2)
        coef_arr = coef_table{w};
        C_table(w) = calc_C(coef_arr, n);
    end
end

function coef_table = make_coef_table(n)
    global nCk
    coef_table = cell(1, floor((n+1)/2));
    for w = 1:floor((n+1)/2)
        sym = gen_monosym(w,n);
        assert(weight(sym, n)==1)
        coef_table{w} = sym;
    end
end

function ent = entropy(coef_arr, n) %having a problem with log of symoblic
    global nCk
    if coef_arr{1} == 0
        ent = sym('0');
    else 
        ent = coef_arr{1}^2 * (2*n-log2(coef_arr{1}^2));
    end
    for k = 1:n
        if coef_arr{k+1} ~= 0
            ent = ent + nCk{n,k} * coef_arr{k+1}^2 * (2*n-log2(coef_arr{k+1}^2));
        end
    end
    ent = ent/(2^(2*n));
end

function inf = influence(coef_arr, n)
    global nCk
    inf = sym('0');
    for k = 1:n
        inf = inf + nCk{n,k} * coef_arr{k+1}^2 * k;
    end
    two = sym('2');
    inf = inf/(two^(2*n));
end

function wt = weight(coef_arr, n)
    global nCk
    wt = coef_arr{1}^2;
    for k = 1:n
        wt = wt + nCk{n, k} * coef_arr{k+1}^2;
    end
    wt = wt/(2^(2*n));
end

function coef_arr = gen_monosym(thresh, n) %WITH IMPLIED 2^-n
    global nCk
    coef_arr = cell(1, n+1);
    num0 = sym('1');
    for m = 1:thresh-1
        num0 = num0 + nCk{n,m};
    end
    two = sym('2');
    coef_arr{1} = 2*num0-two^n;
    for s = 1:n
        numS = krav(thresh - 1, s - 1, n - 1);
        coef_arr{s+1} = 2*numS;
    end
end


function out = krav(k, x, n)
    global nCk
    out = sym(0);
    for j = 0:k
        combo1 = sym(0); %avoids error on n choose k with invalid k
        if j>=0 && j<=x
            if j == 0
                combo1 = sym(1);
            else
                combo1 = nCk{x, j};
            end
        end
        combo2 = sym(0);
        if (k-j)>=0 && (k-j)<=(n-x)
            if k-j == 0
                combo2 = sym(1);
            else
                combo2 = nCk{n-x, k-j};
            end
        end
        prod = combo1 * combo2;
        syms j_sym;
        prod_sign_expr = (-1)^j_sym * prod;
        prod_sign = subs(prod_sign_expr, j_sym, j);
        out = out + prod_sign;
    end
end

function [] = gen_comb_table(n)
    global nCk
    nCk = cell(n,n);
    for m = 1:n
        nCk{m,1} = sym(m);
    end
    for m = 1:n
        for k = 2:m
            if m == k
                nCk{m,k}=1;
            else
                nCk{m,k} = nCk{m-1,k-1} + nCk{m-1, k};
            end
        end
    end
end