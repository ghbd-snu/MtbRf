import lightgbm as lgb
from sklearn.preprocessing import OneHotEncoder
from sklearn.model_selection import train_test_split
from sklearn.metrics import accuracy_score, precision_score, recall_score, roc_auc_score, roc_curve
import pandas as pd

def predict_dst_with_lightgbm(
    variant_data,  # Variant data as a DataFrame
    var_idx_to_use, # List of variant indices to use
    target_data,   # Target data (=resistant to drug) as a Series or array
    split_seed,
    test_size=0.2,
    n_jobs=12,
):
    """
    Train a LightGBM model to predict DST results and evaluate its performance.
    
    Parameters:
    variant_data (pd.DataFrame): DataFrame containing variant data (samples x variants).
    var_idx_to_use (list): List of variant indices to use for the prediction.
    target_data (pd.Series or array): Series or array containing the DST labels.
    split_seed (int): Random seed for splitting the data.
    test_size (float, optional): Proportion of the dataset to include in the test split (default is 0.2).
    n_jobs (int, optional): Number of threads to use for LightGBM (default is 12).
    
    Returns:
    tuple: Accuracy, precision, recall, AUC score, and ROC curve data.
    """
    # Select variants
    X = variant_data.iloc[:, var_idx_to_use].values
    y = target_data.values
    
    # One-hot encoding
    X_train, X_test, y_train, y_test = train_test_split(
        X, y,
        test_size=test_size,
        random_state=split_seed,
    )
    encoder = OneHotEncoder(handle_unknown='ignore')
    encoder.fit(X)
    X_train_encoded = encoder.transform(X_train)
    X_test_encoded = encoder.transform(X_test)
    
    # Convert data into LightGBM's Dataset format
    dtrain = lgb.Dataset(X_train_encoded, label=y_train)
    dtest = lgb.Dataset(X_test_encoded, label=y_test, reference=dtrain)

    params = {
        'objective': 'binary',
        'max_depth': 8,
        'num_leaves': 255,
        'learning_rate': 0.15,
        'colsample_bytree': 0.5,
        'metric': 'binary_logloss',
        'n_jobs': n_jobs,
    }
    model = lgb.train(params, dtrain, num_boost_round=200)
    
    # Make predictions
    predictions = model.predict(X_test_encoded)
    
    # Convert predictions to binary labels (0 or 1)
    predictions_binary = [round(prediction) for prediction in predictions]
    
    # Measure accuracy
    acc = accuracy_score(y_test, predictions_binary)
    prec = precision_score(y_test, predictions_binary)
    rec = recall_score(y_test, predictions_binary)
    auc = roc_auc_score(y_test, predictions)
    roc = roc_curve(y_test, predictions)
    
    return acc, prec, rec, auc, roc
