import http.server
import socketserver
import os

PORT = 8080

web_dir = os.path.join(os.path.dirname(__file__), 'build/html/')
os.chdir(web_dir)

Handler = http.server.SimpleHTTPRequestHandler
httpd = socketserver.TCPServer(("", PORT), Handler)
print("serving at: http://localhost:", PORT, sep="")
httpd.serve_forever()