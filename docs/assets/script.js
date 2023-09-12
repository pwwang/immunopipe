(function($) {
    $('.rst-content a[href^="http"]').attr('target', '_blank');

    if (window.location.hash !== '') {
        $(`.rst-content #faq ~ details${window.location.hash}`).attr('open', true);
    }
})(jQuery);
