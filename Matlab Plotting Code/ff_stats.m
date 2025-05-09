clear 

N = 4;
max_fold = 30;
X=(1:max_fold);

% Each read matrix has 7 columns
fold = zeros(max_fold,7,N);
facets = zeros(30, N);
mileage = zeros(30,N);
% Store the slope of the fitted line
facet_fit_coeffs = zeros(1,N);
for i = 1:N
    file_name = ['f' num2str(i) '.txt'];
    fold(:,:,i) = readmatrix(file_name);
    facets(:,i) = fold(:,2,i);
    mileage(:,i) = fold(:,4,i);
    % Fit line to total facet numbers
    pol = polyfit(X(5:end), facets(5:end,i), 1);
    facet_fit_coeffs(i) = exp(pol(1));
end
%facet_fit_coeffs

% plot facet count stats
for i = 1:N
    if i == 3
        plot(X, facets(:,i),DisplayName=['Fold', num2str(i)],LineWidth=3.5);
    else
        plot(X, facets(:,i),DisplayName=['Fold', num2str(i)]);
    end
    hold on;
end
xlabel("Folds");
ylabel("log(Total Facets)");
legend(Location="northwest");
hold off; 

% plot crease length stats
for i = 1:N
    plot(X, mileage(:,i),DisplayName=['Fold', num2str(i)]);
    hold on;
end
xlabel("Folds");
ylabel("log(Crease Mileage)");
legend(Location="northwest");
hold off; 
