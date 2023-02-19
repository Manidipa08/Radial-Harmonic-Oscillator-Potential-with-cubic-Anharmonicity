clc
clear
n=1001//input("enter number of datapoints : ")
h=6.626D-34
e=1.6D-19 
h_bar=1.054*(1e-34)
l=0 //p wave
b=input("Enter the value of b : ")//J/m^3
b=b*e*(1D+45)*(1D+6)
bb=6
rmin=0//input("enter r_min : ")
rmax=bb*1D-15//input("enter r_max : ")
counter =4
m=1.67*10^(-27)
E0 = 9*1D9
r=linspace(rmin,rmax,n)//range
rnew =r(2:n-1)'
dr=r(2)-r(1)
C=-((h/(2*%pi))^2)/(2*m*(dr^2)*e)
k=100*1.6*(1e-19)*(1e+36)
w=sqrt(k/m)
/*
//rtp calculation
for i=1:4
    rtp(i)=sqrt((2*i-1)*(h/(2*%pi))/(m*w))
end*/
V=zeros(1,n)
for i=1:n
    V(i)=(0.5*k*(r(i)^2)+((1/3)*b*r(i)^3))/e //in eV
end
A=eye(n-2,n-2)
v=diag(V(2:n-1))
//disp(A)
//----------------------By inbuilt command-------------------------
//tic()
D=(-2*C)*ones(n-2,1)
A1=diag(D)
//disp(A1)
C1=C*ones(n-3,1)
A2=diag(C1,1)
//disp(A2)
A3=diag(C1,-1)
//disp(A3)
//Hermitian matrix (tridiagonal)
H=A1+A2+A3+v
//disp("By inbuilt : ",H)
[a1 a2]=spec(H)
Z=spec(H)
for i=1:5
    Eth(i)=(((2*(i-1))+(3/2))*(h/(2*%pi))*w)/e
end
//--------------------------Normalizations & Expectation values-------------------------
vnew = V(2:n-1)
for i = 1:4
    eigenvec = a1(:,i).*a1(:,i)
    Nrm(i)=inttrap(rnew,eigenvec)
    u_r(:,i)=(1/sqrt(Nrm(i))).*a1(:,i)
    neigen(:,i) = u_r(:,i).*u_r(:,i)
    neigen1(:,i) = rnew.*neigen(:,i)
    neigen2(:,i) = rnew.*neigen1(:,i)
    neigen3(:,i) = vnew'.*neigen(:,i)
end
for i=1:counter
    ex_r(i)=inttrap(rnew,neigen1(:,i))
    ex_r2(i)= inttrap(rnew,neigen2(:,i))
    ex_V(i) = inttrap(rnew,neigen3(:,i))
end
for i = 1:n-3
    for j = 1:counter
        r2(i) = (rnew(i)+rnew(i+1))/2
        mid_u(i,j) = (u_r(i,j)+u_r(i+1,j))/2
        diff_u(i,j) = (u_r(i+1,j)-u_r(i,j))/dr
    end
end

for i = 1:n-4
    for j = 1:counter
        r3(i) = (r2(i) + r2(i+1))/2
        mid2_u(i,j) = (u_r(i+2,j) + 2*u_r(i+1,j)+ u_r(i,j))/4
        diff2_u(i,j) = (u_r(i+2,j) - 2*u_r(i+1,j) + u_r(i,j))/(dr*dr)
    end
end
//----------------Uncertainty check--------------------------------------
Un = (4*%pi)/h

for i=1:counter
    y2(:,i) = mid_u(:,i).*diff_u(:,i)
    y3(:,i) = mid2_u(:,i).*diff2_u(:,i)
    ex_p(i) = -1*%i*(h/(2*%pi))*inttrap(r2,y2(:,i))
    ex_p2(i) = -1*(h/(2*%pi))**2*inttrap(r3,y3(:,i))
    sig_r(i) = sqrt(ex_r2(i) - (ex_r(i)*ex_r(i)))
    sig_p(i) = sqrt(ex_p2(i) - (ex_p(i)*ex_p(i)))
    ex_K(i) = ex_p2(i)/(2*m*e)
    un(i) = Un*(sig_r(i).* sig_p(i))
end
//-----------------total energy------------------------
E = ex_V + ex_K
disp("Expectation value <r> =",ex_r(1))
disp("Expectation value <r2> =",ex_r2(1))
disp("Expectation value <p> =",ex_p(1))
disp("Expectation value <p2> =",ex_p2(1))
disp("Expectation value <V> =",ex_V(1))
disp("Expectation value <KE> =",ex_K(1))
disp("Standard deviation of r =",sig_r(1))
disp("Standard deviation of p =",sig_p(1))
disp("Uncertanity Product (hbar/2) = ",un)
disp(" E_comp(eV) from spec         E_comp(eV) (<KE>+<V>) Energy")
disp(string(Z(1))+"        "+string(E(1)))
