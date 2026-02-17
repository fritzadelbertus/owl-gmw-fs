run_alteq:
	python3 alteq.py

test:
	pytest tests/

env:
	python3 -m venv myenv

install:
	pip install -r requirements.txt

clean:
	rm -rf __pycache__/ venv/

check:
	python3 test.py

freeze:
	pip freeze > requirements.txt