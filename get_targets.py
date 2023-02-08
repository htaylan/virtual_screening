from chembl_webresource_client.new_client import new_client

class TargetGetter:
  def get_targets(target_name):
      """
      :param target_name: The name of the inhibitor target
      :return: This function gets targets from ChemBL database
      and return as a dataframe
      """
      self.target = new_client.target
      self.target_query = target.search('target_name')
      self.targets = pd.DataFrame.from_dict(target_query)
      return targets

  def select_target(n):
      """
      :return: Selects the target
      """
      self.selected_target = targets.target_chembl_id[0]
      return selected_target

  def get_inhibitors(selected_target):
      """
      :param selected_target: ChemBl Id of selected target
      :return: The inhibitors for selected target
      from ChemBl database as a dataframe
      """
      self.activity = new_client.activity
      res = activity.filter(target_chembl_id=selected_target)
      df = pd.DataFrame.from_dict(res)
      return df
