# Script for generating params running all 10 folds
import optuna
import pandas as pd
import xgboost as xgb
from sklearn.metrics import mean_absolute_error
import time
import os

# Record the start time
start_time = time.time()

# Change depending on what to run:
file_name = "wing_180k_params.csv"
print("Starting for model ", file_name, "\n")
data = pd.read_feather("C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Data/Processed/wingBV.feather")
df = pd.read_csv('C:/Users/gard_/Documents/MasterThesis/ProjectThesis/MyPipeline/Data/CVfolds/cv_folds_wing_180k.csv')
# Store SNPs in X --- IID 70k, FID 180k
X = data.drop([
            "ID",
            "ringnr",
            "mean_pheno",
            "FID",
            "MAT",
            "PAT",
            "SEX",
            "PHENOTYPE",
            "hatchisland"
        ], axis = 1)
# Some of the SNPS have NA-values. Set to 0
X = X.fillna(0)
# Change from float to int64
X = X.T.astype("int64").T
y = data["ID"]

# X is ringnrs + all SNPs
X_CV = data.drop([
            "ID",
            "mean_pheno",
            "FID",
            "MAT",
            "PAT",
            "SEX",
            "PHENOTYPE",
            "hatchisland"
        ], axis = 1)

# Some of the SNPS have NA-values. Set to 0
X_CV = X_CV.fillna(0)
# Change from float to int64 for all columns not 'ringnr' (i.e. all SNPs)
X_temp = X_CV.drop(['ringnr'], axis = 1)
X_temp = X_temp.T.astype('int64').T
X_temp.insert(0, 'ringnr', X_CV['ringnr'])
X_CV = X_temp

# y is ringnrs + pseudo phenotype
y_CV = data[['ID', 'ringnr']]

for i in range(9,11):
    print("Starting run", i, "at", time.time() - start_time)
    if (i == 10):
        j = 1
    else:
        j = i + 1

    test_idx = df[(df['Fold'] == i) & (df['Set'] == 'test')]['ringnr'].values
    val_idx = df[(df['Fold'] == j) & (df['Set'] == 'test')]['ringnr'].values

    # Training set equal to the intersect of training fold 1 and 2
    f1 = df[((df["Fold"] == i) & (df["Set"] == 'train'))]["ringnr"]
    f2 = df[((df["Fold"] == j) & (df["Set"] == 'train'))]["ringnr"]
    intersect = set(f1).intersection(set(f2))
    # Define training sets
    X_train = X_CV[X_CV["ringnr"].isin(intersect)].drop(["ringnr",], axis = 1)
    y_train = y_CV[y_CV["ringnr"].isin(intersect)].drop(["ringnr",], axis = 1)


    y_val = y_CV[y_CV['ringnr'].isin(val_idx)]["ID"]
    X_val = X_CV[X_CV["ringnr"].isin(val_idx)].drop(["ringnr",], axis = 1)

    def objective(trial):
        """Objective function to be optimized by package optuna.
          Loss-function: MAE.
           n_jobs: Amount of processors used. Computer-specific, needs to be 
                    changed according to computer."""
        params = {
            "objective": "reg:absoluteerror",
            "verbosity": 0,
            "n_estimators": 600, #trial.suggest_int("n_estimators", 50, 300),
            "learning_rate": trial.suggest_float("learning_rate", 1e-3, 0.1, log=True),
            "max_depth": trial.suggest_int("max_depth", 4, 14),
            "subsample": trial.suggest_float("subsample", 0.05, 1.0),
            "colsample_bytree": trial.suggest_float("colsample_bytree", 0.05, 1.0),
            "min_child_weight": trial.suggest_int("min_child_weight", 5, 25),
        }
        model = xgb.XGBRegressor(**params)
        model.fit(X_train, y_train, verbose=False)
        predictions = model.predict(X_val)
        mae = mean_absolute_error(y_val, predictions)
        return mae

    study_full = optuna.create_study(direction='minimize')
    study_full.optimize(objective, n_trials=20)
    df2 = pd.DataFrame([study_full.best_params])
    print("Storing file:\n")
    if os.path.exists(file_name):
        print("Adding to file:\n")
        with open(file_name, mode='a') as f:
            df2.to_csv(f, header=False, index=False)
            f.flush()  # Force the buffer to flush to disk
    else:
        print("Creating file:\n")
        with open(file_name, mode='w') as f:
            df2.to_csv(f, header=True, index=False)
            f.flush()  # Force the buffer to flush to disk

print("Complete!")


