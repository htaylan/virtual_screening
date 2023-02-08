import pandas as pd
from chembl_webresource_client.new_client import new_client


class ChemBLClient:
    def __init__(self):
        self.target = new_client.target
        self.activity = new_client.activity
        self.targets = None
        self.selected_target = None

    def get_targets(self, target_name):
        """
        :param target_name: The name of the inhibitor target
        :return: This function gets targets from ChemBL database
        and return as a dataframe
        """
        target_query = self.target.search(target_name)
        self.targets = pd.DataFrame.from_dict(target_query)
        return self.targets

    def select_target(self, n):
        """
        :return: Selects the target
        """
        self.selected_target = self.targets.target_chembl_id[n]
        return self.selected_target

    def get_inhibitors(self, selected_target=None):
        """
        :param selected_target: ChemBl Id of selected target
        :return: The inhibitors for selected target
        from ChemBl database as a dataframe
        """
        if not selected_target:
            selected_target = self.selected_target
        res = self.activity.filter(target_chembl_id=selected_target)
        df = pd.DataFrame.from_dict(res)
        return df

