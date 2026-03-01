package Cat;
                                     sub new {
                                        my $class = shift;
                                        my $self = {};
                                        bless $self;
                                        if (defined $_[0]) {
                                                $self->{'name'} = shift;
                                                }
                                        if (defined $_[0]) {
                                                $self->{'color'} = shift;
                                                }
                                        return $self;
                                      }

                                     sub meow {
                                        my $self = shift;
                                        print "meow\n";
                                      }

                                     sub printDetails {
                                        my $self = shift;
                                        print "$self->{name}\n";
                                        print "$self->{color}\n";
                                     }

                                     1;
