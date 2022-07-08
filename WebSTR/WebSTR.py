#!/usr/bin/env python3
"""
WebSTR v2 database application
"""

import argparse

#################### Database paths ###############

#################### Set up flask server ###############
server = Flask(__name__)
server.secret_key = 'WebSTR' 

#################### Render region page ###############

# TODO

#################### Render locus page ###############

# TODO

#################### Render HTML pages ###############
@server.route('/')
@server.route('/WebSTR')
def WebSTRHome():
    return render_template('homepage.html')

#################### Set up and run the server ###############
def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument("--host", help="Host to run app", type=str, default="0.0.0.0")
    parser.add_argument("--port", help="Port to run app", type=int, default=5000)
    args = parser.parse_args()
    server.run(debug=False, host=args.host, port=args.port)

if __name__ == '__main__':
    main()