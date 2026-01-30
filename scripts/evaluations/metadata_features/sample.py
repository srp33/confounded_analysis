import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.multioutput import MultiOutputClassifier
from sklearn.inspection import permutation_importance
from sklearn.model_selection import train_test_split
from sklearn.datasets import make_multilabel_classification

# 1. Example data
X, y = make_multilabel_classification(n_samples=1000, n_features=100, n_classes=3, n_labels=6, random_state=42)
feature_names = [f"Feature {i}" for i in range(X.shape[1])]
X_train, X_test, y_train, y_test = train_test_split(X, y, random_state=42)

# 2. Train Model
# MultiOutputClassifier fits one regressor per target label
model = MultiOutputClassifier(HistGradientBoostingClassifier(random_state=42)).fit(X_train, y_train)

# 3. Calculate Per-Label Importance (ROC AUC)
importances = {}
for i, estimator in enumerate(model.estimators_):
    r = permutation_importance(estimator, X_test, y_test[:, i], n_repeats=10, n_jobs = -1, random_state=42, scoring='roc_auc')
    importances[f"Label_{i}"] = r.importances_mean

# 4. Filter & Plot
df_imp = pd.DataFrame(importances, index=feature_names)

# Filter: Keep feature if importance > 0.005 for AT LEAST one label
# You'll store all the data in a csv; filtering can be done later.
df_filtered = df_imp[(df_imp > 0.005).any(axis=1)]

plt.figure(figsize=(10, 18))
sns.heatmap(df_filtered, annot=True, cmap="viridis", fmt=".3f")
plt.title("Feature Importance per Label (ROC AUC)")
plt.show()