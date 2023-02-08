from chembl_webresource_client.new_client import new_client

class getChemBl:
    def __init__(self):
        self.new_client = Client(config.CHEMBL_API_KEY)

    def get_targets(self, target_name):
        """
        :param target_name: The name of the inhibitor target
        :return: This function gets targets from ChemBL database
        and return as a dataframe
        """
        target = self.new_client.target
        target_query = target.search(target_name)
        targets = pd.DataFrame.from_dict(target_query)
        return targets

    def select_target(self, n):
        """
        :return: Selects the target
        """
        selected_target = self.targets.target_chembl_id[n]
        return selected_target

    def get_inhibitors(self, selected_target):
        """
        :param selected_target: ChemBl Id of selected target
        :return: The inhibitors for selected target
        from ChemBl database as a dataframe
        """
        activity = self.new_client.activity
        res = activity.filter(target_chembl_id=selected_target)
        df = pd.DataFrame.from_dict(res)
        return df
