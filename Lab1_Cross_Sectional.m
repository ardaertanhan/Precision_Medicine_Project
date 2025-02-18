%-------------------------------------------------------------------------%
% LAB 1                                                                   %
% Cross-sectional analysis - part 1                                       %
%                                                                         %
% Biomarkers, Precision Medicine & Drug Development                       %
% Prof. Mattia Veronese                                                   %
% Tutors: Marco Pinamonti, Mario Severino                                 %
% Credits: Maria Colpo                                                    %
%-------------------------------------------------------------------------%

clear all
close all
clc

%% 1) Managing the datasets
%% 1.1) Data loading
% Define the path where the .csv files are located
patients_path = fullfile("data", "radiomics_patients.csv");
controls_path = fullfile("data", "radiomics_controls.csv");

% load .csv files as tables
patients = readtable(patients_path);
controls = readtable(controls_path);

%% 1.2) Extract variable names
variable_names_patients = patients.Properties.VariableNames;
variable_names_controls = controls.Properties.VariableNames;

%% 1.3) Tables harmonization
% delete patients row 6 because we are not interested about the Cohort information
patients(:,6) = [];

% change the order of the labels
controls = [controls(:, 1), controls(:, 4:5),controls(:, 2:3),controls(:, 6:end)];

% change name of labels
patients.Properties.VariableNames = controls.Properties.VariableNames;

% change 'Gender' variable to categorical
controls.Gender = categorical(controls.Gender);
patients.Gender = categorical(patients.Gender);


%% 2) Understanding the dataset
%% 2.1) Check demographics matching
% Plots

% Boxplot
figure
boxplot([controls.Age, patients.Age], 'Labels', {'Controls', 'Patients'});
title('Age')
ylabel('Years')

% Barplot
n_gender_c = countcats(controls.Gender);
n_gender_p = countcats(patients.Gender);
X = categorical({'Female', 'Male'});

figure
bar(X, [n_gender_c, n_gender_p])
legend('Controls', 'Patients')
ylabel('#subj')

% Mean & SD
m_age_patients = mean(patients.Age);
m_age_contols = mean(controls.Age);

sd_age_patients = std(patients.Age);
sd_age_contols = std(controls.Age);

disp("%----------------------------------------%")
disp("DEMOGRAPHICS INFO")
disp(['Mean age of patients: ', num2str(m_age_patients)]);
disp(['Mean age of controls: ', num2str(m_age_contols)]);
disp(['SD age of patients: ', num2str(sd_age_patients)]);
disp(['SD age of controls: ', num2str(sd_age_contols)]);
disp('Mean and SD are quite similar so patients and controls are matched')
disp("")

% Pearson Chi-square test
[tbl_match_age, chi2_match_age, p_match_age] = crosstab(patients.Age, controls.Age);
[tbl_match_gender, chi2_match_gender, p_match_gender] = crosstab(patients.Gender, controls.Gender);

disp("%----------------------------------------%")
disp('Pearson Chi-square test:')
disp(['Age match, pval : ', num2str(p_match_age)]);
disp(['Gender match, pval : ', num2str(p_match_gender)]);
disp('both pval>0.05 so patients and controls are matched')
disp("")

%close all

%% 2.2) Data exploration
% Extract numerical arrays from tables
patients_num = table2array(patients(:, [2, 5, 7:end]));
controls_num = table2array(controls(:, [2, 5, 7:end]));
labels_num =  controls(:, [2, 5, 7:end]).Properties.VariableNames;

% Histograms
for ii = 1:size(patients_num, 2)
    figure
    histogram(controls_num(:, ii))
    hold on
    histogram(patients_num(:, ii))
    title(labels_num(ii), 'Interpreter', 'none')
    legend('Controls', 'Patients')
    %pause
end

% Boxplot
for ii = 1:size(patients_num, 2)
    figure
    boxplot([controls_num(:, ii), patients_num(:, ii)], 'Labels', {'Controls', 'Patients'});   
    title(labels_num(ii), 'Interpreter', 'none')
    %pause
end

% Barplot
X = categorical({'Controls', 'Patients'});
for ii = 1:size(patients_num, 2)
    figure
    subplot(211)
    bar(controls_num(:, ii))
    title('Controls')
    subplot(212)
    bar(patients_num(:, ii))
    title('Patients')
    sgtitle(labels_num(ii), 'Interpreter', 'none')
    %pause
end

%close all

%% 2.3) Gaussianity check

% qqplot
for ii = 1:size(patients_num, 2)
    figure
    sgtitle([labels_num(ii)], 'Interpreter', 'none')
    subplot(121)
    qqplot(controls_num(:, ii))
    title('Controls')
    subplot(122)
    qqplot(patients_num(:, ii))
    title('Patients')
    % pause
end

% Lilliefors test
sz = [16, 2];
varTypes = ["int8", "int8"];
varNames = ["controls", "patients"];
rowNames = labels_num;
h_gauss = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames, 'RowNames', rowNames);
p_gauss = table('Size', sz, 'VariableTypes', varTypes, 'VariableNames', varNames, 'RowNames', rowNames);

for ii = 1:size(patients_num, 2)
    [h_gauss.controls(ii), p_gauss.controls(ii)] = lillietest(controls_num(:, ii));
    [h_gauss.patients(ii), p_gauss.patients(ii)] = lillietest(patients_num(:, ii));
end

%close all

%% 2.4) Explore effects of covariates

p_age=[];
chi2_gender = [];
corr_age = [];
p_gender = [];
kk = 1;

for jj = 6:size(controls, 2)
    [corr_age(kk), p_age(kk)] = corr(controls.(jj),controls.Age, 'Type', 'Spearman'); % age
    [tbl_gender, chi2_gender(kk), p_gender(kk)] = crosstab(controls.(jj), controls.Gender); % gender
    
    kk = kk+1;
end

corr_age_table = array2table(corr_age');
label = [controls(:, 6:end).Properties.VariableNames]';
corr_age_table.Properties.RowNames = label;
corr_age_table.Properties.VariableNames = {'Corr Var with Age'};


pval_gender_table = array2table(p_gender');
pval_gender_table.Properties.RowNames = label;
pval_gender_table.Properties.VariableNames = {'Pval Var with Gender'};

for jj = 6:size(controls, 2)
    
    Y_and_X_control = [controls(:, jj), controls(:, 2), controls(:, 3)];
    Y_and_X_control.Gender = categorical(Y_and_X_control.Gender);
    
    var_int_label_control = Y_and_X_control(:, 1).Properties.VariableNames;
    
    fit = fitlm(Y_and_X_control, char(append(var_int_label_control, '~Age*Gender')));
    w = linspace(min(table2array(controls(:, 2))), max(table2array(controls(:, 2))));
    
    figure
    gscatter(table2array(controls(:, 2)), table2array(controls(:, jj)), table2array(controls(:, 3)), [0.8500 0.3250 0.0980; 0.4660 0.6740 0.1880])
    hold on
    line(w, feval(fit, w, 'female'), 'LineWidth', 2, 'Color', [0.8500 0.3250 0.0980])
    line(w, feval(fit, w, 'male'), 'LineWidth', 2, 'Color', [0.4660 0.6740 0.1880])
    title(controls(:, jj).Properties.VariableNames, 'Interpreter', 'none')
    legend('female', 'male')
    xlabel([controls(:, 2).Properties.VariableNames], 'Interpreter', 'none')
    
    %pause
end
%close all

%% 3) Cross Sectional analysis w/o covariates

% 3.1) Already done in point 2.3

% 3.2) Select most appropriate test

% 3.3) Statistical tests comparison

p_vec = zeros(15, 2);

% Wilcoxon-Mann-Whitney
disp('%----------------------%')
disp('Wilcoxon-Mann-Whitney test')
for kk = 6:size(patients, 2)
    
    [p_vec(kk-5, 1),h] = ranksum(table2array(patients(:,kk)), table2array(controls(:,kk)));
    if h == 0
        disp(['Var ', char(patients(:,kk).Properties.VariableNames), ': H0 accepted'])
    else
        disp(['Var ', char(patients(:,kk).Properties.VariableNames), ': H0 rejected'])
    end
end
disp('----------------------')

% t-test
disp('t-test')
for kk = 6:size(patients, 2)
    
    [h, p_vec(kk-5, 2)] = ttest2(table2array(patients(:, kk)), table2array(controls(:, kk)));
    if h == 0
        disp(['Var ', char(patients(:, kk).Properties.VariableNames), ': H0 accepted'])
    else
        disp(['Var ', char(patients(:, kk).Properties.VariableNames), ': H0 rejected'])
    end
end
disp('----------------------')


%% 4) Cross Sectional analysis with covariates

% 4.1) Create new column type
col_to_add = array2table([ones(20,1); 2*ones(20,1)]);
col_to_add.Properties.VariableNames = {'type'};

% 4.2) Build new table controls & patients
new_tab = [controls(:, 1:end); patients(:, 1:end)];
new_tab = [new_tab, col_to_add];
new_tab.type = categorical(new_tab.type);

p_anova = [];
tbl_all_var = struct('Var', {}, 'Pval', {}, 'Table', {});
for jj = 5:(size(new_tab,2)-1)
    
    tbl_all_var(jj).Var = new_tab(:, jj).Properties.VariableNames;
    [tbl_all_var(jj).Pval,tbl_all_var(jj).Table] = anovan(new_tab.(jj), {new_tab.Gender, new_tab.type, new_tab.Age}, 'varnames', {'G' 'T' 'A'}, 'continuous', 3, 'model', [1 0 0; 0 1 0; 0 0 1; 1 0 1]);
    disp(char(tbl_all_var(jj).Var))
    pause
end

close all

