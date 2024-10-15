import os

def get_password():
    password = os.environ.get('TRUBA_PWD')
    if not password:
        raise ValueError("TRUBA_PWD environment variable not set")
    return password


HOSTNAME = "172.16.7.1"
USERNAME = "mtasbas"
HOME_FOLDER = "/truba/home/mtasbas/"
PASSWORD = get_password()
