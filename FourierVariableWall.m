function [F]=FourierVarringWall(y)
%... a function to obtain varriation of via-wall of siw bandpass filter

%% inputs: ================================================================================================

    fcl=12e9;                       %... minimum fequency(GHz) of bandpass(lower cut-off)
    fch=16e9;                       %... maximum fequency(GHz) of bandpass(higher cut-off)
    fmin=10e9;                      %... minimum fequency(GHz) of first half of bandreject
    fmax=18e9;                      %... maximum fequency(GHz) of second half of bandreject
    Df=0.05e9;                      %... delta_f
    Mf=1+(fmax-fmin)/Df;            %... number of fequency samples
    L=20e-3;                        %... length of SIW(milimeter)
    M=80;                           %... number of segments
    N=numel(y);                     %... number of fourier-serie components
    d=0.5e-3;                       %... diameter of vias
    s=0.9e-3;                       %... interspace between vias in each row(center-to-center)
    epsr=3.55;                      %... relative permitivity
    h=0.3048e-3;                    %... substrate thickness(milimeter) 
    kp=1.265;                       %... k-prim constant coefficient in Z0(m) formula
    Zms=50;                         %... matching section impedance(ohm) 
    alpha=30;

%% define variables: ======================================================================================

    Wr=3e8/(2*fcl*sqrt(epsr));
    Dx=L/M;                         %... delta_x
    w=zeros(M,N);
    Weff=zeros(1,M);                %... effective width of SIW
    Wsiw=zeros(1,M);                %... width of SIW
    Zw=zeros(M,Mf);
    Z0=zeros(M,Mf);                 %... charachteristic impedance
    Bet=zeros(M,Mf);                %... beta(phase constant) at SIW
    X=zeros(1,M);
    a=zeros(M,Mf);
    c=zeros(M,Mf);
    ABCD=eye(2,2);                  %... initial value total ABCD matrix(unique)
    E=zeros(1,Mf);                  %... error function
    A=zeros(1,Mf);
    B=zeros(1,Mf);
    C=zeros(1,Mf);
    D=zeros(1,Mf);
    S11=zeros(1,Mf);
    S21=zeros(1,Mf);

%% calculations:    =======================================================================================

    for i=1:Mf
        f=fmin+(i-1)*Df;
        for m=1:M
            X(m)=(m-0.5)*Dx;                                    %... center location of m_th segment along x-axis
            for n=1:N
                w(m,n)=y(n)*cos(2*n*pi*X(m)/L); 
            end 
            Weff(m)=Wr*exp(sum(w(m,:)));                        %... y(n)is Fourier coefficient
            Wsiw(m)=Weff(m)+(d^2)/(0.95*s);
            Bet(i,m)=sqrt((2*pi.*f).^2*epsr*((pi/Weff(m))^2));  
            Zw(i,m)=2*pi.*f./Bet(i,m);
            Z0(i,m)=kp*h.*Zw(i,m)./Weff(m);
            a(m,i)=cos(Dx.*Bet(i));
            c(m,i)=1j.*Z0(m).*sin(Dx.*Bet(m));
            abcd=[a(i,m) c(i,m);c(i,m) a(i,m)];                 %... ABCD matrix of m_th segment(not that here a=d and b=c)
            ABCD=ABCD*abcd;                                     %... Total ABCD matrix
            
        end
        A(i)= ABCD(1,1);
        B(i)= ABCD(1,2); 
        C(i)= ABCD(2,1);
        D(i)= ABCD(2,2);
        S11(i)=(A(i).*Zms+B(i)-C(i).*(Zms^2)-D(i).*Zms)./(A(i).*Zms+B(i)-C(i).*(Zms^2)+D(i).*Zms);
        S21(i)=2*(A(i).*D(i)-B(i).*C(i))/(A(i).*Zms+B(i)-C(i).*(Zms^2)+D(i).*Zms);
        if i>=1+(fcl-fmin)/Df && i<=1+(fch-fmin)/Df
            E(i)=sqrt(alpha.*(abs(S11(i))).^2+(abs(S21(i))-1).^2);
        else
            E(i)=sqrt(alpha.*(abs(S21(i))).^2+(abs(S11(i))-1).^2);
        end 
    end
    
%% object function ========================================================================================

    F=sqrt(sum(E)/Mf);                                          %... objective function
    
end

