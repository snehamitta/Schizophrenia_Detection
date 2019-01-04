labels_train = readtable('train_labels.csv','Delimiter',',');

labels_train.Class = categorical(labels_train.Class, {'Healthy Control','Schizophrenic Patient'}, {0, 1});
summary(labels_train.Class)

FNC_train = dataset('file','train_FNC.csv','Delimiter',',');

summary(FNC_train(:,2:end))
grpstats(cat(2,FNC_train(:,2:end),dataset({labels_train.Class,'Class'})),'Class',{'mean', 'std'})

FNC_test = dataset('file','test_FNC.csv','Delimiter',',');


SBM_train = dataset('file','train_SBM.csv','Delimiter',',');

summary(SBM_train(:,2:end))

grpstats(cat(2,SBM_train(:,2:end),dataset({labels_train.Class,'Class'})),'Class',{'mean', 'std'})
SBM_test = dataset('file','test_SBM.csv','Delimiter',',');

example = dataset('file','submission_example.csv','Delimiter',',');

export(example,'file','new_submission.csv','Delimiter',',');