%% 等额本息的月供计算公式
function monthlyPayment = calculateEqualMonthlyPayment(principal, annualInterestRate, loanYears)
    % 计算月利率
    monthlyInterestRate = annualInterestRate / 12 / 100;
    
    % 计算总的还款月数
    totalMonths = loanYears * 12;
    
    % 使用等额本息公式计算月还款金额
    monthlyPayment = principal * monthlyInterestRate * (1 + monthlyInterestRate)^totalMonths / ((1 + monthlyInterestRate)^totalMonths - 1);
end
