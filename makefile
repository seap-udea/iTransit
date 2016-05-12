DATABASE=Aristarchus
BACKDIR=data/dump

clean:
	find . -name "*~" -exec rm {} \;

cleandata:
	@rm -r data/Aristarco6/*

cleanall:clean cleandata

commit:
	@echo "Commiting changes..."
	@-git commit -am "Commit"
	@git push origin master

pull:
	@echo "Pulling from repository..."
	@git reset --hard HEAD	
	@git pull
	@chown -R www-data.www-data .

backup:
	@echo "Backuping sinfin..."
	@bash backup.sh 

restore:
	@echo "Restoring table $TABLE..."
	@-p7zip -d $(BACKDIR)/$(DATABASE).tar.7z
	@-tar xf $(BACKDIR)/$(DATABASE).tar
	@echo -n "Enter root mysql password: "
	@mysql -u root -p $(DATABASE) < $(BACKDIR)/$(DATABASE).sql
	@p7zip $(BACKDIR)/$(DATABASE).tar

permissions:
	@echo "Setting web permissions..."
	@chown -R www-data.www-data .

edit:
	@emacs -nw makefile *.php web/*.php css/*.css js/*.js
