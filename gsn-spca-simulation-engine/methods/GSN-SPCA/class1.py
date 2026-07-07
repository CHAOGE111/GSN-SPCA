import pandas as pd

# Read the second column of a.csv
a_df = pd.read_csv('a.csv')
second_column = a_df.iloc[:, 1]

# Read the file GSE174330_norm_counts_FPKM_GRCh38.p13_Symbol.csv
gse_df = pd.read_csv('GSE174330_norm_counts_FPKM_GRCh38.p13_Symbol.csv')

# Extract the rows corresponding to the values in the second column of a.csv
rows_to_extract = gse_df.iloc[second_column - 1]  # Subtract 1 because Python indexing starts at 0

# Create a new DataFrame, the first column is the second column of a.csv, and the rest are the corresponding row data
new_df = pd.concat([second_column.reset_index(drop=True), rows_to_extract.reset_index(drop=True)], axis=1)

# Output to a new file, separated by spaces
new_df.to_csv('output.txt', sep=' ', index=False, header=False)


import pandas as pd

# Assume the data file is 'output.txt' and separated by spaces
df = pd.read_csv("output.txt", sep=" ", header=None)

# Transpose the data
df_transposed = df.T

# Reset index and set column names for the transposed DataFrame
df_transposed.reset_index(inplace=True)
df_transposed.columns = ['Name'] + [f'Feature_{i+1}' for i in range(df_transposed.shape[1] - 1)]

# Prepare the target variable column with initial values as None or another placeholder
df_transposed['Target'] = None

# Assume rows 2-51 (index 2–51) belong to class 0, rows 52-101 belong to class 1
df_transposed.loc[2:51, 'Target'] = 0  
df_transposed.loc[52:101, 'Target'] = 1  

# Save the transformed data to a new CSV file
df_transposed.to_csv('transformed_data.csv', index=False)


# 3
import pandas as pd
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score
from sklearn.metrics import confusion_matrix
from sklearn.metrics import precision_score
from sklearn.metrics import recall_score
from sklearn.metrics import f1_score
from sklearn.metrics import roc_auc_score

# Load data
data = pd.read_csv('transformed_data.csv')

# Drop the second row
data = data.drop(labels=data.index[1], axis=0)

# Separate features and target
X = data.iloc[:, 1:-1]  # Features are from 2nd column to the second-to-last column
y = data.iloc[:, -1]    # Target is the last column

# Check and handle NaN values in the target variable
if y.isnull().values.any():
    # Drop rows with NaN in the target
    data = data.dropna(subset=[y.name])
    # Update X and y
    X = data.iloc[:, 1:-1]
    y = data.iloc[:, -1]

# Split train and test sets
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Create Random Forest Classifier
rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)

# Train the model
rf_classifier.fit(X_train, y_train)

# Predict on test set
y_pred = rf_classifier.predict(X_test)

# Calculate accuracy
accuracy = accuracy_score(y_test, y_pred)

# Get feature importance
importances = rf_classifier.feature_importances_

# Save feature importance in a DataFrame
feature_importance_df = pd.DataFrame({'Feature': X.columns, 'Importance': importances})

# Export feature importance table to CSV file without sorting
feature_importance_df.to_csv('feature_importance.csv', index=False)

# Calculate other evaluation metrics
conf_matrix = confusion_matrix(y_test, y_pred)
precision = precision_score(y_test, y_pred, average='macro')
recall = recall_score(y_test, y_pred, average='macro')
f1 = f1_score(y_test, y_pred, average='macro')


# 4
# Read the first column of output.txt
output_columns = pd.read_csv('output.txt', sep='\s+', header=None)[0]

# Load feature_importance.csv
feature_importances = pd.read_csv('feature_importance.csv')

# Reset index and insert 'gene' column at the first position
feature_importances.reset_index(drop=True, inplace=True)  
feature_importances.insert(0, 'gene', output_columns)  

# Assume N = 1.5
N = 1.5

# Normalize importance values
if 'Importance' in feature_importances.columns:
    min_importance = feature_importances['Importance'].min()
    max_importance = feature_importances['Importance'].max()
    feature_importances['Importance_normalized'] = 1 + (N - 1) * (feature_importances['Importance'] - min_importance) / (max_importance - min_importance)

# Save updated feature importance to a new CSV file
feature_importances.to_csv('updated_feature_importance.csv', index=False)


import pandas as pd

# Read updated_feature_importance.csv (comma separated) and result-1_p2.txt (tab separated)
feature_importance_df = pd.read_csv('updated_feature_importance.csv', delimiter=',')
result_df = pd.read_csv('result-1_p2.txt', delimiter='\t', header=None)

# Use the first column as index
feature_importance_df.set_index(feature_importance_df.columns[0], inplace=True)

# Define accuracy threshold
Q = 0.6  
# Function to replace values based on Importance_normalized
def replace_value(row):
    new_row = []
    for value in row:
        # If value exists in the index, replace with Importance_normalized
        if value in feature_importance_df.index:
            new_row.append(feature_importance_df.at[value, 'Importance_normalized'])
        else:
            new_row.append(1)
    return pd.Series(new_row)

# Apply the replacement function to each row
result_df = result_df.apply(replace_value, axis=1)

# Take the maximum value of each row
result_df['max_value'] = result_df.max(axis=1)

# If accuracy < Q, set all max_value to 1
if accuracy < Q:
    result_df['max_value'] = 1

# Save the results to a txt file
result_df['max_value'].to_csv('result_max_values.txt', index=False, header=False)
print("ok")
