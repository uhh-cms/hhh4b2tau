from order import Config
from columnflow.config_util import add_category 

def add_categories(
    cfg: Config,
):

    add_category(
            cfg,
            id=1,
            name="incl",
            selection="cat_incl",
            label="inclusive",
        )
    add_category(
        cfg,
        name="2j",
        selection="cat_2j",
        label="2 jets",
    )
