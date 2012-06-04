#!/usr/local/bin/perl -w

################################################################################
# N.B. Run "perl configuration-assistant.perl" if you see a message like:
# configuration-assistant.pl: /bin/perl: bad interpreter: No such file or directory
################################################################################
my $VERSION = '1.0.5';
################################################################################

use strict;

use Cwd "getcwd";
use Fcntl ':mode';
use File::Basename qw(basename dirname);
use File::Path "mkpath";
use File::Spec;
use IO::Handle;
use Getopt::Long "GetOptions";

sub println { print @_; print "\n"; }

STDOUT->autoflush(1);

my $R = 1;
my $W = 2;
my $X = 4;
my $CREATED = 8;
my $FOUND = 16;
my $RWX = $FOUND + $R + $W + $X;

println "==========================================";
println "The SRA toolkit documentation page is";
println "http://www.ncbi.nlm.nih.gov/Traces/sra/std";
println "==========================================\n";

my %options;
Help() unless (GetOptions(\%options, 'help', 'version', 'wordy'));
Help() if ($options{help});
Version() if ($options{version});

$options{references} = 1 if (defined @ARGV);

my ($VDB_CONFIG, $BIN_DIR) = find_bin("vdb-config");

my @KFG_FILES;
{
    print "checking configuration files... ";
    my $f = `$VDB_CONFIG -f`;
    if ($f =~ /\w/) {
        while (chomp $f) {}
        @KFG_FILES = split /\n/, $f;
        println "'$f'";
    } else {
        println "not found";
    }
}

my %refseq_cfg;
%refseq_cfg = check_refseq_kfg() if (@KFG_FILES);

umask 0077;

my $HOME_CFG_DIR;
my $fixed_config = do_refseq_cfg();
$fixed_config   |= do_krypto_cfg();
do_schema_cfg();

if ((! defined $options{references}) && $fixed_config == 0 && 0)
{   $options{references} = 1 }
$options{references} = askyN
       ("Would you like to test cSRA files for remote reference dependencies?",
        $fixed_config)
    unless ($options{references});
exit 0 unless ($options{references});

my $CHECK_ONLY = $refseq_cfg{refseq_dir} eq 'READ^ONLY';

my ($ALIGN_INFO, $ALIGN_INFO_DIR) = find_bin("align-info");
WARN("using align-info and vdb-config located in different directories:\n"
        . "\t$ALIGN_INFO\n\t$VDB_CONFIG")
    if ($BIN_DIR ne $ALIGN_INFO_DIR);

my $WGET;
$WGET = find_wget() unless ($CHECK_ONLY);

my $got_packed;
my @packed;

if ($#ARGV > -1) {
    foreach (@ARGV) {
        load($CHECK_ONLY, $refseq_cfg{refseq_dir}, $_);
    }
} else {
    while (1) {
        my $f = ask("Enter cSRA file name (Press Enter to exit)");
        last unless ($f);
        load($CHECK_ONLY, $refseq_cfg{refseq_dir}, $f);
    }
}

################################################################################
# FUNCTIONS #
################################################################################

sub do_refseq_cfg {
    my $fixed_config = 0;

    if ($refseq_cfg{NO_volumes} || $refseq_cfg{NO_servers}) {
        fix_refseq_cfg(\%refseq_cfg,
            "Existing configuration is invalid. Whould you like to fix it?",
            "Your tools will not work properly.\n" .
            "For more information visit SRA website.");
    }

    if ($refseq_cfg{refseq_dir}) {
        if ($refseq_cfg{FIX_volumes} || $refseq_cfg{FIX_paths}) {
            $refseq_cfg{refseq_dir}
                = ask_refseq_change($refseq_cfg{refseq_dir});
            if ($refseq_cfg{refseq_dir} ne 'READ^ONLY') {
                my $config_dir
                    = del_all_refseq_cfg_path($refseq_cfg{refseq_dir});
                create_refseq_cfg
                    ($config_dir, $refseq_cfg{refseq_dir}, $BIN_DIR);
            }
        } elsif (not $refseq_cfg{refseq_dir_prm}) {
            unless (askYn("Configuration is found\n"
                . "but the local reference repository $refseq_cfg{refseq_dir} "
                . "does not exist.\n"
                . "Do you want to create it?"))
            {
                println "Your tools will not work properly.";
                println "Remote reference sequences will not be resolved.";
                println "For more information visit SRA website.";
                exit 1;
            }
            my $prm = check_dir($refseq_cfg{refseq_dir}, "create_missing");
            my $created = $prm & $CREATED;
            $prm &= ~ $CREATED;
            if ($prm != $RWX) {
                if ($^O ne 'cygwin' || $created != $CREATED) {
                    FATAL("Cannot create $refseq_cfg{refseq_dir}");
                    fix_refseq_cfg(\%refseq_cfg,
                        "Local reference repository is read-only. ".
                        "Whould you like to change it?",
                        "You will not be able to upload reference sequences.\n"
                        . "For more information visit SRA website.");
                }
            }
            ++$fixed_config;
        }
    } else {
        unless (askYn("Do you want to create a new configuration?")) {
            println "Your tools will not work properly without configuration.";
            println "For more information visit SRA website.";
            exit 1;
        }

        print "checking configuration directories... ";
        my $paths = `$VDB_CONFIG -d`;
        if ($?) {
            die $!;
        }
        while (chomp $paths) {}
        my @paths = split /:/, $paths;
        println if ($#paths >= 0);

        my $config_dir;
        foreach (@paths) {
            next if ($_ eq '/etc/ncbi');
            if ($^O eq 'MSWin32') { # Windows: translate POSIX to Windows path
                tr|/|\\|;
                s/^\\([a-zA-Z])\\/$1:\\/;
            } elsif ($^O eq 'cygwin') {
                if (m|^/([a-zA-Z])/|) {
                    $_ = "/cygdrive$_";
                }
            }

            my $prm = check_dir($_, "create_missing");
            my $created = $prm & $CREATED;
            $prm &= ~ $CREATED;
            if ($prm == $RWX || ($^O eq 'cygwin' and $created == $CREATED)) {
                $config_dir = $_;
                last;
            }
        }
        FATAL("cannot find a usable directory to create configuration file")
            unless ($config_dir);

        ($refseq_cfg{refseq_dir}, $refseq_cfg{refseq_dir_prm})
            = make_refseq_dir();

        create_refseq_cfg($config_dir, $refseq_cfg{refseq_dir}, $BIN_DIR);
        ++$fixed_config;
    }

    if ($refseq_cfg{refseq_dir} ne 'READ^ONLY') {
        clean_empty($refseq_cfg{refseq_dir});
    }

    return $fixed_config;
}

sub ask {
    my ($prompt) = @_;
    print "$prompt: ";
    my $in = <STDIN>;
    unless ($in) {
        println;
        return "";
    }
    chomp $in;
    return $in;
}
sub askYn { return askYN($_[0], 'yes');}
sub askyN { return askYN($_[0], 'no' );}
sub askYN {
    my ($q, $yes) = @_;
    $yes = '' if ($yes eq 'no');
    print "$q ";
    if ($yes) {
        print "[Y/n] ";
    } else {
        print "[y/N] ";
    }
    my $in = <STDIN>;
    chomp $in;
    if ($in) {
        return $in eq 'Y' || $in eq 'y'
          || $in eq 'YES' || $in eq 'yes' || $in eq 'Yes';
    } else {
        return $yes;
    }
}
sub ask_refseq_change {
    my ($refseq) = @_;
    die unless ($refseq);
    my $force;
    while (1) {
        println "Your repository directory is $refseq.";
        println "It is read-only. You cannot add new sequences to it.";
        println "Make your choice:";
        println "1) Use existing repository";
        println "2) Use a different repository";
        println "3) exit the script";
        print "Your selection? [1] ";
        my $in = <STDIN>;
        chomp $in;
        return 'READ^ONLY' if (! $in || $in eq "1");
        if ($in eq "2") {
            last;
        }
        exit 0 if ($in eq "3");
    }
    my ($path, $perm) = make_refseq_dir();
    return $path;
}
sub make_refseq_dir {
    my $deflt;
    if ($^O eq 'cygwin' and $ENV{USERPROFILE}) {
        $deflt = $ENV{USERPROFILE};
        $deflt =~ tr|\\|/|;
        $deflt =~ s|^([a-zA-Z]):/|/$1/|;
        if ($1) {
            $deflt = "/cygdrive$deflt";
        }
    } elsif ($ENV{HOME}) {
        $deflt = $ENV{HOME};
    } elsif ($ENV{USERPROFILE}) {
        $deflt = $ENV{USERPROFILE};
    } else {
        $deflt = getcwd();
    }
    $deflt = "." unless($deflt);
    $deflt = File::Spec->catdir($deflt, "ncbi", "refseq");

    while (1) {
        my $path;
        print "Specify installation directory for reference objects";
        if ($deflt) {
            print " [$deflt]";
        }
        print ": ";

        my $in = <STDIN>;
        unless ($in) {
            println;
            $path = "";
        }
        chomp $in;
        if ($in) {
            $path = $in;
        } elsif ($deflt) {
            $path = $deflt;
        }
        exit 1 unless ($path);
        my $perm = check_dir($path, "create_missing");
        my $prm = $perm;
        my $created = $prm & $CREATED;
        $prm &= ~ $CREATED;
        if ($prm == $RWX || ($^O eq 'cygwin' and $created == $CREATED)) {
            return ($path, $perm);
        } elsif ($prm & $FOUND) {
            println "'$path' does not seem to be writable.";
            println "You will not be able to add new sequences to it.";
            if (askyN("Do you want to use it?")) {
                return ($path, $perm);
            }
        }
    }
}

sub POSIXify {
    ($_) = @_;
    # convert to POSIX path
    s|^/cygdrive/|/|;
    tr|\\|/|;
    s|^([a-zA-Z]):/|/$1/|;
    return $_;
}

sub create_refseq_cfg {
    my ($config_dir, $refseq_dir, $bin_dir) = @_;

    $config_dir = mk_home_cfg_dir();

    my $kfg;
    for (my $i = 0; ; ++$i) {
        $kfg = "$config_dir/ncbi-config";
        $kfg .= "$i" if ($i);
        $kfg .= ".kfg";
        last unless (-e $kfg);
    }

    print "creating configuration file $kfg... ";
    $refseq_dir = POSIXify($refseq_dir);

    open(F, ">$kfg") or die "cannot open $kfg";
    print F "refseq/paths = \"$refseq_dir\"\n";
    close F or die "cannot close $kfg";
    println "ok";

    my %tmp_cfg = check_refseq_kfg();
    die unless ($tmp_cfg{refseq_dir_prm} = $RWX);
}

sub Version {
      my $bin = basename($0);
      print << "END";
$bin version $VERSION
END

      exit 1;
}
sub Help {
      my $bin = basename($0);
      print << "END";
$bin version $VERSION

Utility to help configure the SRA tools to be able
to access the local reference repository,
determine which reference sequences a cSRA file relies
upon and to fetch them from NCBI.
Fetched references are placed in the local reference repository.

Usage: $bin [--wordy] [cSRA_FILE...]

Options
-v, --version     print version and exit
-w, --wordy       increase "fetcher" verbosity
END

      exit 1;
}

sub fix_refseq_cfg {
    my ($refseq_cfg, $question, $bye) = @_;
    unless (askYn($question)) {
        print $bye;
        exit 1;
    }
    my $config_dir = del_all_refseq_cfg_path($refseq_cfg);
    ($refseq_cfg->{refseq_dir}, $refseq_cfg{refseq_dir_prm})
        = make_refseq_dir();
    create_refseq_cfg($config_dir, $refseq_cfg->{refseq_dir}, $BIN_DIR);
}

sub del_all_refseq_cfg_path {
    my ($refseq_cfg) = @_;
    my $config_dir;
    my %bak;
    foreach (@KFG_FILES) {
        my $f = $_;
        my $file_bak;
        my $fixed = "";
        open IN, $f or die "cannot open $_";
        while (<IN>) {
            if (m{refseq/servers} || m{refseq/volumes} || m{refseq/paths}) {
                $config_dir = dirname($f) unless ($config_dir);
                $fixed .= "# $_";
                for (my $i = 0; not $file_bak; ++$i) {
                    $file_bak = "$f.bak";
                    $file_bak .= "$i" if ($i);
                    last unless (-e $file_bak);
                }
            } else {
                $fixed .= $_;
            }
        }
        close $f;
        if ($file_bak) {
            unless (rename($f, $file_bak)) {
                rollback(%bak);
                FATAL("Cannot update(save) $f: is file/directory writable?");
            } else {
                $bak{$file_bak} = $f;
                unless (open(OUT, ">$f")) {
                    rollback(%bak);
                    FATAL("Cannot update $f: $!");
                }
                print OUT $fixed;
                unless (close(OUT)) {
                    rollback(%bak);
                    FATAL("Cannot update(close) $f: $!");
                }
            }
        }
    }
    return $config_dir;
}

sub rollback {
    my (%bak) = @_;
    while (my ($bak, $orig) = each %bak) {
        rename($bak, $orig);
    }
}

sub do_schema_cfg {
    my $error;
    print "checking schema configuration... ";
    my $tmp = `$VDB_CONFIG vdb/schema/paths 2>&1`;
    if ($? == 0) {
        chomp $tmp;
        if ($tmp =~ m|<paths>(.+)</paths>|) {
            my $paths = $1;
            println $paths;
            my @paths = split(/:/, $paths);
            $error += !check_schema_file('align/align.vschema', @paths);
            $error += !check_schema_file('ncbi/seq.vschema'   , @paths);
            $error += !check_schema_file('vdb/vdb.vschema'    , @paths);
        } else {
            $error = "unexpected: $tmp";
            println $error;
        }
    } elsif ($tmp =~ /path not found/) {
        $error = "not found";
        println $error;
    } else {
        $error = "unknown vdb-config schema error";
        println $error;
    }
    if ($error) {
        println "--------------------------------------";
        println "WARNING: SCHEMA FILES CANNOT BE FOUND.";
        println "IT COULD CAUSE LOADERS TO FAIL.";
        println "--------------------------------------";
    }
    return ! $error;
}

sub check_schema_file {
    my ($file, @dir) = @_;
    print "checking $file... ";
    foreach my $dir(@dir) {
        my $path = "$dir/$file";
        if (-e $path) {
            println $path;
            return 1;
        }
    }
    println "not found";
    return 0;
}

sub load {
    my ($check_only, $refseq_dir, $f) = @_;
    println "Determining $f external dependencies...";
    my $cmd = "$ALIGN_INFO $f";
    my @info = `$cmd`;
    if ($?) {
        println "$f: error";
    } else {
        my $refs = 0;
        my $ok = 0;
        my $ko = 0;
        foreach (@info) {
            chomp;
            my @r = split /,/;
            if ($#r >= 3) {
                my ($seqId, $remote) = ($r[0], $r[3]);
                my $refseqName;
                ++$refs;
                if ($remote eq 'remote' && ! $check_only) {
                    ($refseqName, $remote) = check_packed($seqId, $remote);
                }
                if ($remote eq 'remote') {
                    if ($check_only) {
                        println "$seqId: unresolved";
                        ++$ok;
                    } else {
                        my $res = download($refseqName);
                        if ($res) {
                            ++$ko;
                        } else {
                            ++$ok;
                        }
                    }
                }
            }
        }
        if ($check_only) {
            if ($ok == 0) {
                println "All $refs references were found";
            } else {
                println "$refs references were checked, $ok missed";
            }
        } else {
            print "All $refs references were checked (";
            print "$ko failed, " if ($ko);
            println "$ok downloaded)";
        }
    }
}

sub download {
    my ($file) = @_;
    my $refseq_dir = $refseq_cfg{refseq_dir};
    print "Downloading $file... ";
    println if ($options{wordy});
    my $dest = "$refseq_dir/$file";
    my $cmd = "$WGET \"$dest\"" .
        " http://ftp-private.ncbi.nlm.nih.gov/sra/refseq/$file 2>&1";
    my $res;
    if ($options{wordy}) {
        $res = system($cmd);
    } else {
        `$cmd`;
        $res = $?;
    }
    print "$file: " if ($options{wordy});
    if ($res) {
        println "failed";
        if (-z $dest)
        {   unlink($dest); }
    } else {
        println "ok";
    }
    return $res;
}

sub check_packed {
    my ($seqId, $remote) = @_;
    my $file = $seqId;
    my $refseq_dir = $refseq_cfg{refseq_dir};
    unless ($got_packed) {
        download('packed.txt');
        ++$got_packed;
        my $f = "$refseq_dir/packed.txt";
        if (-e $f) {
            if (open(IN, $f)) {
                foreach (<IN>) {
                    next if (/^\s*\#/);
                    chomp;
                    if (/^\s*([\w\.]+)\s*$/) {
                        push (@packed, $1);
                    }
                }
                close(IN);
            }
        }
    }
    foreach (@packed) {
        if ($seqId =~ /^$_/) {
            $file = $_;
            if (-s "$refseq_dir/$file") {
                $remote = 'local';
            }
            last;
        }
    }
    return ($file, $remote);
}

sub clean_empty {
    my ($refseq_dir) = @_;
    print "checking $refseq_dir for invalid empty reference files... ";
    my $i = 0;
    opendir DIR, $refseq_dir or die "cannot opendir $refseq_dir";
    while ($_ = readdir DIR) {
        next if (/^\.{1,2}$/);
        my $f = "$refseq_dir/$_";
        my $empty;
        if (-z $f) {
            ++$empty;
        } elsif (-s $f < 999) {
            open F, $f or die "cannot open $f";
            my $data = '';
            while (<F>) {
                while(chomp) {};
                $data .= $_;
            }
            ++$empty if ($data =~ m|<title>404 - Not Found</title>|);
        }
        if ($empty) {
            unlink $f or die "cannot remove $f";
            ++$i;
        }
    }
    closedir DIR;
    if ($i)
    {   println "$i found"; }
    else
    {   println "none found"; }
}

################################################################################

sub check_dir {
    my ($dir, $create_missing) = @_;

    $dir = File::Spec->canonpath($dir);
    print "checking $dir... ";

    my $prm = 0;
    unless (-e $dir) {
        println "not found";
        return 0 unless ($create_missing);

        print "checking ${dir}'s parents... ";
        $dir = File::Spec->canonpath($dir);
        my @dirs = File::Spec->splitdir($dir);
        my $test = File::Spec->rootdir();
        if ($^O eq 'MSWin32') {
            $test = "";
        } else {
            FATAL("bad root directory '$test'") unless (-e $test);
        }
        foreach (@dirs) {
            my $prev = $test;
            if ($test) {
                $test = File::Spec->catdir($test, $_);
            } else {
                $test = File::Spec->catdir($_, File::Spec->rootdir());
            }
            if (! -e $test) {
                $test = $prev;
                last;
            }
        }

        print "($test)... ";
        my $cygwin_beauty;
# cygwin does not detect -r for $ENV{USERPROFILE}
        if (! -r $test || ! -x $test) {
            if ($^O eq 'cygwin') {
                ++$cygwin_beauty;
            } else {
                println "not readable";
                return 0;
            }
        }
        if (! -x $test) {
            if ($^O eq 'cygwin') {
                ++$cygwin_beauty;
            } else {
                println "not writable";
                return 0;
            }
        }
        if ($cygwin_beauty) {
            println("fail to check");
        } else {
            println("ok");
        }

        print "creating $dir... ";
        unless (mkpath($dir)) {
            die "cannot mkdir $dir" unless ($cygwin_beauty);
            println "failed. Is it writable?";
            return 0;
        }
        println("ok");
        $prm += $CREATED;
        print "checking $dir... ";
    }

    $prm += $FOUND;

    {
        my $cygwin_beauty;
        my $failures;
        if (-r $dir) {
            $prm += $R;
        }
        if (-w $dir) {
            $prm += $W;
        }
        if (-x $dir) {
            $prm += $X;
        }

        if (! -r $dir || ! -x $dir) {
            if ($^O eq 'cygwin') {
                ++$cygwin_beauty;
            } else {
                println "not readable";
                ++$failures;
            }
        }
        if (! $failures and ! -w $dir) {
            if ($^O eq 'cygwin') {
                ++$cygwin_beauty;
            } else {
                println "not writable";
                ++$failures;
            }
        }
        if ($cygwin_beauty) {
            println("fail to check");
        } elsif (!$failures) {
            println("ok");
        }
    }

    return $prm;
}

sub check_refseq_kfg {
    print "checking refseq configuration... ";
    my %refseq_cfg;
    refseq_from_kfg(\%refseq_cfg);

    if ($refseq_cfg{servers} and !$refseq_cfg{volumes}) {
        ++$refseq_cfg{NO_volumes};
        println "found invalid";
        return %refseq_cfg;
        FATAL("need to update existing invalid configuration");
        # TODO
    } elsif ($refseq_cfg{volumes} and !$refseq_cfg{servers}) {
        ++$refseq_cfg{NO_servers};
        println "found invalid";
        FATAL("need to update existing invalid configuration");
        # TODO
    } elsif ($refseq_cfg{servers} and $refseq_cfg{volumes}) {
        println "servers=$refseq_cfg{servers} volumes=$refseq_cfg{volumes}";
        if ($refseq_cfg{paths}) {
           WARN("/refseq/servers&volumes found: /refseq/paths will be ignored");
        }
    } elsif ($refseq_cfg{paths}) {
        println "paths=$refseq_cfg{paths}";
    } else {
        println "not found";
        return %refseq_cfg;
    }

    if ($refseq_cfg{servers} and index($refseq_cfg{servers}, ":") != -1) {
        die "Unexpected configuration: servers=$refseq_cfg{servers}";
    }
    if ($refseq_cfg{volumes} and index($refseq_cfg{volumes}, ":") != -1) {
        die "Unexpected configuration: volumes=$refseq_cfg{volumes}";
    }
    if ($refseq_cfg{paths} and index($refseq_cfg{paths}, ":") != -1) {
        die "Unexpected configuration: paths=$refseq_cfg{paths}";
    }

    my $dir;
    if ($refseq_cfg{servers} and $refseq_cfg{volumes}) {
        $dir = File::Spec->catdir($refseq_cfg{servers}, $refseq_cfg{volumes});
    } else {
        $dir ="$refseq_cfg{paths}";
    }

    if ($^O eq 'MSWin32') { # Windows: translate POSIX to Windows path
        $dir =~ tr|/|\\|;
        $dir =~ s/^\\([a-zA-Z])\\/$1:\\/;
    } elsif ($^O eq 'cygwin' and $dir =~ m|^/[a-zA-Z]/|) {
        $dir = "/cygdrive$dir";
    }

    $refseq_cfg{refseq_dir} = $dir;
    my $prm = check_dir($dir);
    $refseq_cfg{refseq_dir_prm} = $prm;
    if ($prm == 0) { # not found
        return %refseq_cfg;
    } elsif ($prm != $RWX) {
        if ($^O ne 'cygwin') {
            if ($refseq_cfg{servers}) {
                ++$refseq_cfg{FIX_volumes};
            } else {
                ++$refseq_cfg{FIX_paths};
            }
#           FATAL("refseq repository is invalid or read-only");
        } # else cygwin does not always can tell permissions correctly
    }
    return %refseq_cfg;
}

sub find_wget {
    my $WGET;

    print "checking for wget... ";
    my $out = `wget -h 2>&1`;
    if ($? == 0) {
        println "yes";
        $WGET = "wget -O";
    } else {
        println "no";
    }

    unless ($WGET) {
        print "checking for curl... ";
        my $out = `curl -h 2>&1`;
        if ($? == 0) {
            println "yes";
            $WGET = "curl -o";
        } else {
            println "no";
        }
    }
    unless ($WGET) {
        print "checking for ./wget... ";
        my $cmd = dirname($0) ."/wget";
        my $out = `$cmd -h 2>&1`;
        if ($? == 0) {
            println "yes";
            $WGET = "$cmd -O";
        } else {
            println "no";
        }
    }
    FATAL("none of wget, curl could be found") unless ($WGET);

    return $WGET;
}

sub find_bin {
    my ($name) = @_;

    my $prompt = "checking for $name";
    my $basedir = dirname($0);

    # built from sources
    print "$prompt (local build)... ";
    if (-e File::Spec->catfile($basedir, "Makefile")) {
        my $f = File::Spec->catfile($basedir, "build", "Makefile.env");
        if (-e $f) {
            my $dir = `make -s bindir -C $basedir 2>&1`;
            if ($? == 0) {
                chomp $dir;
                my $try = File::Spec->catfile($dir, $name);
                print "($try";
                if (-e $try) {
                    print ": found)... ";
                    my $tmp = `$try -h 2>&1`;
                    if ($? == 0) {
                        println "yes";
                        return ($try, $dir);
                    } else {
                        println "\nfailed to run '$try -h'";
                    }
                } else {
                    println ": not found)";
                }
            }
        }
    } else {
        println "no";
    }

    # try the script directory
    {
        my $try = File::Spec->catfile($basedir, $name);
        print "$prompt ($try";
        if (-e $try) {
            print ": found)... ";
            my $tmp = `$try -h 2>&1`;
            if ($? == 0) {
                println "yes";
                return ($try, $basedir);
            } else {
                println "\nfailed to run '$try -h'";
            }
        } else {
            println ": not found)";
        }
    }

    # check from PATH
    {
        my $try = "$name";
        print "$prompt ($try)... ";
        my $tmp = `$try -h 2>&1`;
        if ($? == 0) {
            println "yes";
            return ($try, "");
        } else {
            println "no";
        }
    }

    FATAL("$name could not be found");
}

sub do_krypto_cfg {
    my $fixed_config = 0;
    print "checking krypto configuration... ";
    my $nm = 'pwfile';
    my $name = "krypto/$nm";
    my $v = `$VDB_CONFIG $name 2>&1`;
    if ($?) {
        die $! unless ($v =~ /path not found while opening node/);
        println "not found";
        if (askyN("Do you want to set VDB password?")) {
            update_krypto_cfg('create kfg');
            $fixed_config = 1;
        }
    } else {
        $v =~ /<$nm>(.*)<\/$nm>/;
        die "Invalid '$nm' configuration" unless ($1);
        $v = $1;
        println "$nm=$v";
        if (askyN("Do you want to update VDB password?")) {
            update_krypto_cfg();
            $fixed_config = 1;
        }
    }
    return $fixed_config;
}

sub mk_home_cfg_dir {
    return $HOME_CFG_DIR if ($HOME_CFG_DIR);
    my $config_dir;
    if ($ENV{HOME}) {
        $config_dir = $ENV{HOME};
    } elsif ($ENV{USERPROFILE}) {
        $config_dir = $ENV{USERPROFILE};
    } else { die 'none of $HOME, $USERPROFILE is defined' }
    print "checking local configuration directory... ";
    $config_dir = File::Spec->catdir($config_dir, '.ncbi');
    print "$config_dir... ";
    if (-e $config_dir) {
        println "ok";
        return $HOME_CFG_DIR = $config_dir;
    } else {
        println "not found";
        print "creating $config_dir... ";
        die "cannot mkdir $config_dir" unless (mkpath($config_dir));
        println("ok");
        return $HOME_CFG_DIR = $config_dir;
    }
}

sub update_krypto_cfg {
    my ($update_kfg) = @_;
    my $VDB_PASSWD = (find_bin('vdb-passwd'))[0];
    if ($update_kfg) {
        my $config_dir = mk_home_cfg_dir();
        my $pwfile = File::Spec->catdir($config_dir, '.vdbpass');
        $pwfile = POSIXify($pwfile);
        my $kfg;
        for (my $i = 0; ; ++$i) {
            $kfg = "$config_dir/krypto";
            $kfg .= "$i" if ($i);
            $kfg .= ".kfg";
            last unless (-e $kfg);
        }
        print "creating configuration file $kfg... ";
        open(F, ">$kfg") or die "cannot open $kfg";
        print F "krypto/pwfile = \"$pwfile\"\n";
        close F or die "cannot close $kfg";
        println "ok";
    }
    my $res = system("$VDB_PASSWD -q");
    if ($res) {
        println "password update failed";
    } else {
        println "password updated ok";
    }
}

sub refseq_config {
    my ($nm) = @_;
    my $v = `$VDB_CONFIG refseq/$nm 2>&1`;
    if ($?) {
        if ($v =~ /path not found while opening node/) {
            $v = '';
        } else {
            die $!;
        }
    } else {
        $v =~ /<$nm>(.*)<\/$nm>/;
        die "Invalid 'refseq/$nm' configuration" unless ($1);
        # TODO die if (refseq/paths = "") : fix it
        $v = $1;
    }
    return $v;
}

sub refseq_from_kfg {
    my ($refseq) = @_;
    $refseq->{servers} = refseq_config('servers');
    $refseq->{volumes} = refseq_config('volumes');
    $refseq->{paths} = refseq_config('paths');
}

sub FATAL {
    my ($msg) = @_;
    print basename($0);
    println ": FATAL: $msg";
    exit 1;
}

sub WARN {
    my ($msg) = @_;
    print basename($0);
    println ": WARNING: $msg";
}

################################################################################
# EOF #
################################################################################
