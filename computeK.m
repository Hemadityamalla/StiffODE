function kr = computeK(order,eqNo,i,f,h,y,x,J)
    kr = zeros(order,eqNo);
    gamma = 1.0 + 1/sqrt(2.0); %chosen such that the scheme is stable for larger step sizes.
    dx = 1e-5; %Used to estimate the derivative
    b = zeros(eqNo,1);
    dfdt = zeros(eqNo,1);
        %First slope- k1
    for ii=1:eqNo
        dfdt(ii,1) = (f{ii}(x(i-1)+dx,y(i-1,:)) - f{ii}(x(i-1),y(i-1,:)))/dx;%Maybe replace this by the derivative function in DERIVSUITE 
        b(ii,1) = f{ii}(x(i-1),y(i-1,:)); %Extra term arises for a non-autonomous case
    end
    kr(1,:) = (eye(length(J)) - h*gamma*J)\b;
        %Second slope- k2
    for ii=1:eqNo
        b(ii,1) = f{ii}(x(i-1)+h,y(i-1,:)+h*kr(1,ii)) - 2*kr(1,ii); %Extra term arises for the non-autonomous case.
    end
    kr(2,:) = (eye(length(J)) - h*gamma*J)\b;
    
end