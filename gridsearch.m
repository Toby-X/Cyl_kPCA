function [besttheta] = gridsearch(theta1,theta2,input,score,phi,mu,para_test,test)
%GRIDSEARCH for parameters in kriging based on rmse/mean relative error
%   Detailed explanation goes here
errmin = 1e8;
besttheta = [theta1(1),theta2(1)];

for m = 1:length(theta1)
    for n = 1:length(theta2)
        theta = [theta1(m),theta2(n)];
        [dmodel,~] = dacefit(input,score,@regpoly0,@corrgauss,theta);
        
        [beta_predict] = predictor(para_test,dmodel);
        pred=phi*beta_predict'+repmat(mu,1,size(beta_predict,1));
        
        rmse=sqrt(mean((test - pred).^2,'all'));
        disp(rmse)
        
%         mean_rel_err1=0;
%         num_order=10.^(floor(log10(abs(max(val1(:,specie))))));
%         sum=0;
%         for p=1:length(pred1)
%             if val1(p,specie)>0.05*num_order
%                 rel_err=abs(val1(p,specie)-pred1)/(val1(p,specie)+ 0.001*num_order);
%                 mean_rel_err1=mean_rel_err1+rel_err;
%                 sum=sum+1;
%             end
%         end
%         mean_rel_err1=mean_rel_err1*100/sum;  %以百分数表示
% 
%         mean_rel_err2=0;
%         num_order=10.^(floor(log10(abs(max(val2(:,specie))))));
%         sum=0;
%         for p=1:length(pred2)
%             if val2(p,specie)>0.05*num_order
%                 rel_err=abs(val2(p,specie)-pred2)/(val2(p,specie)+ 0.001*num_order);
%                 mean_rel_err2=mean_rel_err2+rel_err;
%                 sum=sum+1;
%             end
%         end
%         mean_rel_err2=mean_rel_err2*100/sum;  %以百分数表示
%         mean_rel_err = mean([mean_rel_err1, mean_rel_err2]);
        
        if rmse < errmin
            errmin = rmse;
            besttheta = theta;
        end
    end
end

end

