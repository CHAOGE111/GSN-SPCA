import pandas as pd

# ------------------------------------------------------------
# Step 1: Extract rows from gse_df based on indices in a.csv
# ------------------------------------------------------------

# Read the second column from a.csv
a_df = pd.read_csv('a.csv')
second_column = a_df.iloc[:, 1]

# Read main dataset
gse_df = pd.read_csv('data.csv')

# Extract rows from gse_df using indices in second_column
# (subtract 1 because Python indexing starts at 0)
rows_to_extract = gse_df.iloc[second_column - 1]

# Create a new dataframe: first column is from a.csv, 
# followed by corresponding rows from gse_df
new_df = pd.concat([second_column.reset_index(drop=True), 
                    rows_to_extract.reset_index(drop=True)], axis=1)

# Save to new file with space as the delimiter
new_df.to_csv('output.txt', sep=' ', index=False, header=False)

# ------------------------------------------------------------
# Step 2: Transpose the extracted data and prepare labels
# ------------------------------------------------------------

df = pd.read_csv("output.txt", sep=" ", header=None)

# Transpose data
df_transposed = df.T

# Reset index and assign column names
df_transposed.reset_index(inplace=True)
df_transposed.columns = ['Name'] + [f'Feature_{i+1}' for i in range(df_transposed.shape[1] - 1)]

# Add target column (initially None)
df_transposed['Target'] = None

# Assign class labels (0 or 1) to specific row ranges
df_transposed.loc[2:97, 'Target'] = 0  
df_transposed.loc[98:195, 'Target'] = 1  

# Save transformed dataset
df_transposed.to_csv('transformed_data.csv', index=False)

# ------------------------------------------------------------
# Step 3: Train a Random Forest classifier and compute importance
# ------------------------------------------------------------

from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, confusion_matrix, precision_score, recall_score, f1_score

# Load transformed data
data = pd.read_csv('transformed_data.csv')

# Drop the second row
data = data.drop(labels=data.index[1], axis=0)

# Separate features and target
X = data.iloc[:, 1:-1]  
y = data.iloc[:, -1]    

# Handle missing values in target
if y.isnull().values.any():
    data = data.dropna(subset=[y.name])
    X = data.iloc[:, 1:-1]
    y = data.iloc[:, -1]

# Train-test split
X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)

# Train Random Forest classifier
rf_classifier = RandomForestClassifier(n_estimators=100, random_state=42)
rf_classifier.fit(X_train, y_train)

# Predictions
y_pred = rf_classifier.predict(X_test)

# Compute accuracy
accuracy = accuracy_score(y_test, y_pred)

# Compute feature importance
importances = rf_classifier.feature_importances_
feature_importance_df = pd.DataFrame({'Feature': X.columns, 'Importance': importances})

# Save feature importance (unsorted)
feature_importance_df.to_csv('feature_importance.csv', index=False)

# Additional evaluation metrics
conf_matrix = confusion_matrix(y_test, y_pred)
precision = precision_score(y_test, y_pred, average='macro')
recall = recall_score(y_test, y_pred, average='macro')
f1 = f1_score(y_test, y_pred, average='macro')

# ------------------------------------------------------------
# Step 4: Normalize feature importance and update results
# ------------------------------------------------------------

# Read gene IDs from output.txt
output_columns = pd.read_csv('output.txt', sep='\s+', header=None)[0]

# Load feature importance
feature_importances = pd.read_csv('feature_importance.csv')
feature_importances.reset_index(drop=True, inplace=True)

# Insert gene column
feature_importances.insert(0, 'gene', output_columns)

# Normalization factor
N = 4  

# Normalize importance values
if 'Importance' in feature_importances.columns:
    min_importance = feature_importances['Importance'].min()
    max_importance = feature_importances['Importance'].max()
    feature_importances['Importance_normalized'] = 1 + (N - 1) * (
        (feature_importances['Importance'] - min_importance) / (max_importance - min_importance)
    )

# Save updated importance file
feature_importances.to_csv('updated_feature_importance.csv', index=False)

# ------------------------------------------------------------
# Step 5: Replace values in result-1_p2.txt with normalized importance
# ------------------------------------------------------------

# Load updated importance and result file
feature_importance_df = pd.read_csv('updated_feature_importance.csv', delimiter=',')
result_df = pd.read_csv('result-1_p2.txt', delimiter='\t', header=None)

# Use gene column as index
feature_importance_df.set_index(feature_importance_df.columns[0], inplace=True)

# Define accuracy threshold
Q = 0.6  

# Function to replace values with normalized importance
def replace_value(row):
    new_row = []
    for value in row:
        if value in feature_importance_df.index:
            new_row.append(feature_importance_df.at[value, 'Importance_normalized'])
        else:
            new_row.append(1)
    return pd.Series(new_row)

# Apply replacement
result_df = result_df.apply(replace_value, axis=1)

# Extract max value per row
result_df['max_value'] = result_df.max(axis=1)

# If accuracy is below threshold, override max_value with 1
if accuracy < Q:
    result_df['max_value'] = 1

# Save results
result_df['max_value'].to_csv('result_max_values.txt', index=False, header=False)
print("ok")
