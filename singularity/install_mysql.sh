### Installing MySQL, credits: https://serverfault.com/a/783528
export DEBIAN_FRONTEND=noninteractive

MYSQL_ROOT_PASSWORD="my_secret_password"

echo debconf mysql-server/root_password password $MYSQL_ROOT_PASSWORD | debconf-set-selections
echo debconf mysql-server/root_password_again password $MYSQL_ROOT_PASSWORD | debconf-set-selections
apt-get -qq install mysql-server=8.0.19-0ubuntu5 > /dev/null

### Install Expect and dependancies
cd /opt
wget http://cz.archive.ubuntu.com/ubuntu/pool/universe/e/expect/tcl-expect_5.45.4-2build1_amd64.deb
wget http://cz.archive.ubuntu.com/ubuntu/pool/universe/e/expect/expect_5.45.4-2build1_amd64.deb
wget http://security.ubuntu.com/ubuntu/pool/main/t/tzdata/tzdata_2022a-0ubuntu0.20.04_all.deb
wget http://cz.archive.ubuntu.com/ubuntu/pool/main/t/tcl8.6/libtcl8.6_8.6.10+dfsg-1_amd64.deb
wget http://cz.archive.ubuntu.com/ubuntu/pool/main/t/tcl8.6/tcl8.6_8.6.10+dfsg-1_amd64.deb

dpkg -i tzdata_2022a-0ubuntu0.20.04_all.deb
dpkg -i libtcl8.6_8.6.10+dfsg-1_amd64.deb
dpkg -i tcl8.6_8.6.10+dfsg-1_amd64.deb
dpkg -i tcl-expect_5.45.4-2build1_amd64.deb
dpkg -i expect_5.45.4-2build1_amd64.deb

    # Build Expect script
tee ~/secure_our_mysql.sh > /dev/null << EOF
spawn $(which mysql_secure_installation)

expect "Enter password for user root:"
send "$MYSQL_ROOT_PASSWORD\r"

expect "Press y|Y for Yes, any other key for No:"
send "y\r"

expect "Please enter 0 = LOW, 1 = MEDIUM and 2 = STRONG:"
send "2\r"

expect "Change the password for root ? ((Press y|Y for Yes, any other key for No) :"
send "n\r"

expect "Remove anonymous users? (Press y|Y for Yes, any other key for No) :"
send "y\r"

expect "Disallow root login remotely? (Press y|Y for Yes, any other key for No) :"
send "y\r"

expect "Remove test database and access to it? (Press y|Y for Yes, any other key for No) :"
send "y\r"

expect "Reload privilege tables now? (Press y|Y for Yes, any other key for No) :"
send "y\r"

EOF

    # Run Expect script, this runs the "mysql_secure_installation" script which removes insecure defaults.
expect ~/secure_our_mysql.sh

    # Cleanup
rm -v ~/secure_our_mysql.sh 

echo "MySQL setup completed. Insecure defaults are gone. Please remove this script manually when you are done with it (or at least remove the MySQL root password that you put inside it."