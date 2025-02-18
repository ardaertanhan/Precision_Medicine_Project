clear all;
clc;

%---------------------1)Managing the datasets--------------------------- 

% Hastalar veri dosyasını bir tabloya yükle
patientsData = readtable('radiomics_patients.csv');

% Kontroller veri dosyasını bir tabloya yükle
controlsData = readtable('radiomics_controls.csv');

% Hastalar tablosunun değişken adlarını yazdır
disp('Patients Variable Names:');
disp(patientsData.Properties.VariableNames);

% Kontroller tablosunun değişken adlarını yazdır
disp('Controls Variable Names:');
disp(controlsData.Properties.VariableNames);

% Her iki tablonun da sütun adlarını al
patientsVariables = patientsData.Properties.VariableNames;
controlsVariables = controlsData.Properties.VariableNames;

% Hastalar tablosunda olup, kontroller tablosunda olmayan değişkenleri bul
patientsUniqueVariables = setdiff(patientsVariables, controlsVariables);

% Kontroller tablosunda olup, hastalar tablosunda olmayan değişkenleri bul
controlsUniqueVariables = setdiff(controlsVariables, patientsVariables);

% Bulunan değişkenleri ekrana yazdır
disp('Patients Unique Variables:');
disp(patientsUniqueVariables);

disp('Controls Unique Variables:');
disp(controlsUniqueVariables);

%---------------------2)Understanding the dataset-----------------------

% 1. Demographics (age and gender) matching for patients vs controls
%    Check if demographics are matched by comparing distributions

% Demographics for patients
patientsAge = patientsData.age;
patientsSex = patientsData.sex;

% Demographics for controls
controlsAge = controlsData.age;
controlsSex = controlsData.sex;

% Plot age distributions for patients and controls
figure;
subplot(1,2,1);
histogram(patientsAge, 'Normalization', 'probability', 'FaceColor', 'blue');
xlabel('Age');
ylabel('Probability');
title('Patients Age Distribution');
subplot(1,2,2);
histogram(controlsAge, 'Normalization', 'probability', 'FaceColor', 'red');
xlabel('Age');
ylabel('Probability');
title('Controls Age Distribution');
legend('Controls');

% Explore gender distribution for patients and controls
figure;
subplot(1,2,1);
pie(categorical(patientsSex));
title('Patients Gender Distribution');
subplot(1,2,2);
pie(categorical(controlsSex));
title('Controls Gender Distribution');

% 2. Explore the data with the most appropriate plot
%    You may need to decide what type of plot is most suitable based on the data characteristics

% For example, if you have numerical variables, you can use histograms or boxplots
% For categorical variables, you can use bar plots or pie charts

% 3. Check the Gaussianity of the distributions
%    You can use statistical tests like Shapiro-Wilk test or visual inspection using histograms

% For example, perform Shapiro-Wilk test for age distributions
[H_patientsAge, p_patientsAge, ~] = swtest(patientsAge);
[H_controlsAge, p_controlsAge, ~] = swtest(controlsAge);
disp(['Shapiro-Wilk Test p-value for Patients Age: ', num2str(p_patientsAge)]);
disp(['Shapiro-Wilk Test p-value for Controls Age: ', num2str(p_controlsAge)]);

% 4. Explore effects of covariates: compare each target variable with demographics data in controls only
%    Use appropriate statistical tests like t-test or ANOVA for numerical variables, and chi-square test for categorical variables

% For example, compare a target variable with age in controls using t-test
% Spearman korelasyonu hesapla
rho_age_gender = corr(patientsData.age, double(categorical(patientsData.sex)), 'Type', 'Spearman');

% Gender'a göre yaş dağılımını gösteren çapraz-tablo oluştur
crosstab_age_gender = crosstab(categorical(patientsData.sex), patientsData.age);

% Korelasyon sonucunu ekrana yazdır
disp(['Spearman Correlation between Age and Gender: ', num2str(rho_age_gender)]);

% Çapraz-tabloyu ekrana yazdır
disp('Cross-tabulation of Age by Gender:');
disp(crosstab_age_gender);

% % Cinsiyet değişkenini kategorik hale getir
% patientsData.sex = categorical(patientsData.sex);
% 
% mdl_male = fitlm(patientsData(patientsData.sex == '0', :).age, patientsData(patientsData.sex == '0', :).morph_diam);
% mdl_female = fitlm(patientsData(patientsData.sex == '1', :).age, patientsData(patientsData.sex == '1', :).morph_diam);
% 
% % gscatter kullanarak grafiği çiz
% figure;
% gscatter(patientsData.sex, patientsData.morph_diam, patientsData.sex);
% hold on;
% 
% % Erkekler için regresyon çizgisini çiz
% x_male = linspace(0.75, 1.25, 100); % Erkeklerin cinsiyet kodu 0, biraz sağa kaydırıyoruz
% y_male = predict(mdl_male, zeros(size(x_male))); % Erkeklerin regresyon çizgisi
% plot(x_male, y_male, 'r-', 'LineWidth', 2);
% 
% % Kadınlar için regresyon çizgisini çiz
% x_female = linspace(1.75, 2.25, 100); % Kadınların cinsiyet kodu 1, biraz sağa kaydırıyoruz
% y_female = predict(mdl_female, ones(size(x_female))); % Kadınların regresyon çizgisi
% plot(x_female, y_female, 'b-', 'LineWidth', 2);
% 
% % Grafik için genel ayarlar
% xlabel('Gender');
% ylabel('Morph Diam');
% title('Regression Lines by Gender');
% legend('Male', 'Female', 'Location', 'best');
% grid on;
% hold off;

%-----------------------------------------------------------LAB2---


% Değişkenlerin adları

for i = 1:length(variableNames)
    % Belirlenen test türü
    testType = determineTestType(patientsData.(variableNames{i}));
    
    % Testi uygula
    switch testType
        case 'ranksum'
            [p_value, ~, stats] = ranksum(patientsData.(variableNames{i}), controlsData.(variableNames{i}));
        case 'ttest2'
            [h, p_value, ~, stats] = ttest2(patientsData.(variableNames{i}), controlsData.(variableNames{i}));
    end
    
    % Sonuçları ekrana yazdır
    disp(['Comparison between Patients and Controls for ', variableNames{i}, ':']);
    disp(['p-value: ', num2str(p_value)]);
    if exist('h', 'var')
        disp(['t-statistic: ', num2str(stats.tstat)]);
    end
end

function testType = determineTestType(data)
    % Determine the appropriate test type based on the characteristics of the data
    
    % Check for normality using Shapiro-Wilk test
    [~, p, ~] = swtest(data);
    
    % Determine test type based on normality and sample size
    if p < 0.05 || length(data) < 30
        % If not normal or sample size is small, use Wilcoxon-Mann-Whitney test
        testType = 'ranksum';
    else
        % Otherwise, use T-test
        testType = 'ttest2';
    end
end


