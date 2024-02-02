from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import classification_report
import pickle
import pandas as pd
from sklearn.base import BaseEstimator
import numpy as np


class BaseModel:
    def __init__(self, estimator:BaseEstimator=RandomForestClassifier()) -> None:
        """
        Initialize the object.

        Parameters
        ----------
        estimator: BaseEstimator
            A scikit-learn estimator
        
        Returns
        -------
            None
        """
        self.estimator = estimator

    def fit(self, X:pd.DataFrame, y:pd.Series) -> None:
        """
        Fit the ML model estimator

        Parameters
        ----------
        X: pd.DataFrame
            Pandas dataframe containing the training sequence aminoacid content.
        Y: pd.Series
            Pandas series containing the training labels
    
        Returns
        -------
            None
        """
        self.estimator.fit(X, y)
    
    def predict(self, X:pd.DataFrame) -> np.array:
        """
        Generates a prediction based on the fitted model.

        Parameters
        ----------
        X: pd.DataFrame
            Protein sequence test dataset 

        Returns
        -------
        y_pred: np.array
            Predicted y value based on the model
        """
        return self.estimator.predict(X)
        
    def validate(self, X_test:pd.DataFrame, y_test:pd.Series) -> str:
        """
        Validates the model using test data

        Parameters
        ----------
        X_test: pd.DataFrame
            Pandas testing dataframe with sequence data

        y_test: pd.series
            Pandas testing dataframe with labels data
        Returns
        -------
        classification_report: str
            Report with the main classification metrics
        """
        y_pred = self.predict(X_test)
        return classification_report(y_test, y_pred)
    
    def save(self, file_path:str) -> None:
        """
        Saves the model to a serialized pickle file.

        Parameters
        ----------
        file_path: str
            path to the model file
        
        Returns
        -------
            None
        """
        with open(file_path, 'wb') as handle:
            pickle.dump(self, handle)

